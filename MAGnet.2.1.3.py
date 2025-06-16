#! /env/cns/proj/agc/scratch/conda/miniconda3/bin/python

import os,sys,re
import argparse
import csv
import json
from collections import defaultdict, Counter
from Bio import SeqIO
import sqlite3
from sqlite3 import Error
from multiprocessing import Pool
from functools import partial
import logging
import time

def create_table(db_file, sql_statements):
    """ create a table from the sql_statements
    :param db_file: SQLite database to connect (creates it if it doesn't exist)
    :param sql_statements: CREATE TABLE statements
    """
    try:
        with sqlite3.connect(db_file) as conn:
            cursor = conn.cursor()
            for statement in sql_statements:
                cursor.execute(statement)
            conn.commit()
    except Error as e:
        print(e)

def data2insert(infos2insert, keys_order):
    """ move key, value of a dictionary into a list of tuples
    :param infos2insert: dictionary with object of the table as key and all informations as values
    :param keys_order: the ordered columns of the table
    :return: list of tuples that include values to insert 
    """
    list_data =[]
    for data, infos in infos2insert.items():
        tuple_infos = [data]
        for field in keys_order:
            if field in infos:
                value = infos[field]
                if value is None:
                    tuple_infos.append("NULL")
                elif isinstance(value, list):
                    tuple_infos.append(",".join(map(str, value)))
                else:
                    tuple_infos.append(value)
            else:
                tuple_infos.append("NULL")
        list_data.append(tuple(tuple_infos))    
    return list_data

def add_data(conn, sql_statement, data):
    """ insert a new row into the a table
    :param conn: Database connection object
    :param sql_statement: INSERT INTO table statement
    :param data: tuple (or list) that include values to insert
    """
    cur = conn.cursor()
    cur.execute(sql_statement, data)
    conn.commit()
    return cur.lastrowid


def populate_table(db_file, sql_statement, list_data):
    """ insert values into a table
    :param db_file: SQLite database to connect
    :param sql_statement: INSERT INTO table statement
    :param list_data: list of tuples that include values to insert 
    """
    try:
        with sqlite3.connect(db_file) as conn:
            for data in list_data:
                data_id = add_data(conn, sql_statement, data)
                #print(f'Created a new record with the id {data_id}')
                
    except sqlite3.Error as e:
        print(e)

def get_contigs_info(db_file, sql_query):
    """  query all rows from the sql_query statement
    :param db_file: ANVIO database (SQLite database to connect)
    :param sql_query: SELECT statement
    :return: list of tuples (all rows returned by the query) or None
    """
    try:
        with sqlite3.connect(db_file) as conn:
            cur = conn.cursor()
            cur.execute(sql_query)  
            contigs_info = cur.fetchall()
        
            return contigs_info
    except sqlite3.Error as e:
        print(e)
        return None
    
def get_bins_info(db_file, sql_query, collection):
    """  query all rows from the sql_query statement
    :param db_file: ANVIO PROFILE.db (SQLite database to connect)
    :param sql_query: SELECT statement
    :param collection: ANVIO collection
    :return: list of tuples (all rows returned by the query) or None
    """
    try:
        with sqlite3.connect(db_file) as conn:
            cur = conn.cursor()
            cur.execute(sql_query, (collection,))  
            profile_info = cur.fetchall()
            #print(profile_info)
            return profile_info
    except sqlite3.Error as e:
        print(e)
        return None

def get_GTDB_taxonomy(gtdb_metadata_filename):
    """ extract the taxonomy of each representative genome in the GTDB  
    :param gtdb_metadata_filename: path to the GTDB metadata for all bacterial and archaea genomes
    :return: dictionary of GTDB representative accession to GTDB taxonomy
    """
    access2taxo = defaultdict(dict)
    #cpt=0
    with open(gtdb_metadata_filename) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)    
        for line in reader :
            access_nb = line[0]
            representative = line[18]
            taxo = line[19]
            if representative == "t":
                 #cpt+=1
                 access2taxo[access_nb]=taxo
            #if cpt==5000:
                #break
    print(f"length access2taxo: {len(access2taxo)}")
    return access2taxo

def merge_dicts(dict1, dict2):
    merged_data = {}
    all_keys = set(dict1.keys()) | set(dict2.keys())
    
    for key in all_keys:
        if key in dict1.keys():
            merged_data[key] = dict1.get(key, {})
        elif key in dict2.keys():
            merged_data[key] = dict2.get(key, {})
    
    return merged_data

def get_GTDB_proteins(gtdb_protein_filename):
    """ give the correspendance between one gene and genome access in the GTDB
    :param gtdb_path: path to the gtdb current release
    :return: dictionary with gene identifier as key and GTDB genome access as value
    """
    orf2access={}
    with open(gtdb_protein_filename,"r") as file:
        for line in file :
            line=line.rstrip()
            if line[0] == ">" :
                liste=line.split(" ")
                access = liste[0]
                access = access[1:]
                orf = liste[1]
                orf2access[orf]=access

    return orf2access

def run_diamond_blastp(diamond_path, gtdb_protein_db, assembly_protein_filename, cpu):
    """ run diamond blastp for all single GTDB genome
    :param diamond_path: path to the diamond directory
    :param gtdb_protein_db: path to the gtdb proteins database 
    :param assembly_protein_filename: path to the assembly proteins 
    :param cpu: number of threads necessary to run diamond
    """
######  warning add module load diamond to the sbatch ############

    diamond_filename = os.path.join(diamond_path, f"diamond.out")
    cmd = "diamond blastp --db " + gtdb_protein_db + " --query " + assembly_protein_filename + " --out " + diamond_filename + " --mid-sensitive --outfmt 6 sseqid qseqid pident evalue bitscore --threads " + str(cpu)

    status = os.system(cmd)

    if status != 0:
        sys.exit(f"[Error] something went wrong with diamond for {access}, exit.")
    else:
        done = wdir + "/done"
        with open(done, 'w') as done_file:
            done_file.write('Diamond completed successfully\n')

def get_diamond_besthits(diamond_path):
    """ give diamond best hits (e_value < 0.001) for all the assembly genes
    :param diamond_path: path to the diamond output directory
    :return: dictionary of assembly gene id to best diamond scores
    """
    asm_orf2best_hits={}
    diamond_filename = os.path.join(diamond_path, f"diamond.out")
    if not os.path.exists(diamond_filename):
        print(f"Warning: File {diamond_filename} not found...")
        return
        
    with open(diamond_filename,"r") as f:
        for line in f:
            line = line.rstrip()
            liste = line.split("\t")
            gtdb_orf=liste[0]
            asm_orf=liste[1]
            e_value=float(liste[3])
            bitscore=float(liste[4])
            if e_value < 0.001 :
                if asm_orf not in asm_orf2best_hits :
                    asm_orf2best_hits[asm_orf] = {}
                if bitscore not in asm_orf2best_hits[asm_orf]:
                    asm_orf2best_hits[asm_orf][bitscore]=[]                         
                asm_orf2best_hits[asm_orf][bitscore].append(gtdb_orf)

    #os.remove(diamond_filename)
    return asm_orf2best_hits

def get_consensus_taxonomy(BESTHITS_LIST) :
    """ give a consensus taxonomy from a list of taxonomy
    :param BESTHITS_LIST: list of taxonomy
    :return: consensus taxonomy as STRING
    """
    taxo_gene=""
    species_set=set()
    genus_set=set()
    family_set=set()
    order_set=set()
    class_set=set()
    phylum_set=set()
    domain_set=set()
    for taxo2 in BESTHITS_LIST :
        liste=taxo2.split(";")
        domain_set.add(liste[0])
        phylum_set.add(liste[1])
        class_set.add(liste[2])
        order_set.add(liste[3])
        family_set.add(liste[4])
        genus_set.add(liste[5])
        species_set.add(liste[6])
    if len(domain_set) ==1 :
        taxo_gene=list(domain_set)[0]
    else :
        taxo_gene="d__"
    if len(phylum_set) ==1 :
        taxo_gene=taxo_gene + ";" + list(phylum_set)[0]
    else :
        taxo_gene=taxo_gene + ";p__"
    if len(class_set) ==1 :
        taxo_gene=taxo_gene + ";" + list(class_set)[0]
    else :
        taxo_gene=taxo_gene + ";c__"
    if len(order_set) ==1 :
        taxo_gene=taxo_gene + ";" + list(order_set)[0]
    else :
        taxo_gene=taxo_gene + ";o__"
    if len(family_set) ==1 :
        taxo_gene=taxo_gene + ";" + list(family_set)[0]
    else :
        taxo_gene=taxo_gene + ";f__"
    if len(genus_set) ==1 :
        taxo_gene=taxo_gene + ";" + list(genus_set)[0]
    else :
        taxo_gene=taxo_gene + ";g__"
    if len(species_set) ==1 :
        taxo_gene=taxo_gene + ";" + list(species_set)[0]
    else :
        taxo_gene=taxo_gene + ";s__"

    return taxo_gene

def k_best_hits(bitscore2gtdb_orf,k,access2taxo):
    """ give the consensus gene taxonomy from the k diamond best hits  
    :param bitscore2gtdb_orf: dictionnary for diamond bitscore to GTDB protein
    :param k: number of best hits to keep
    :param orf2access: dictionnary for GTDB protein to GTDB accession  
    :param access2taxo: dictionnary for GTDB accession to GTDB taxonomy  
    :return: assembly gene consensus taxonomy as STRING
    """
    cpt = 0
    besthit_list=[]
    for bitscore in sorted(bitscore2gtdb_orf,reverse=True):
        for access_orf in bitscore2gtdb_orf[bitscore]:
            cpt+=1
            access,orf=access_orf.split("__")
            taxo = access2taxo[access]

        if cpt <= k :
            besthit_list.append(taxo)
            #print('\t'+str(bitscore) +'\t'+ orf +'\t'+ access +'\t'+taxo)
        else:
            break
 
    return get_consensus_taxonomy(besthit_list)

def get_gene_taxonomy(anvio_contig_list, asm_orf2best_hits, k, access2taxo):
    """ give a list of the genes taxonomy for each contig of the assembly
    :param anvio_contig_list: anvio contigs list
    :param asm_orf2best_hits: dictionary for assembly's protein to best diamond scores
    :param k: number of best hits to keep
    :param access2taxo: dictionary for GTDB accession to GTDB taxonomy
    :return: dictionary for contig to gene taxonomy as LIST
    """
    contig2taxo={}    
    for asm_orf in asm_orf2best_hits:
        #print(asm_orf)
        gene_taxo=k_best_hits(asm_orf2best_hits[asm_orf],k,access2taxo)
        #print('\t'+gene_taxo)
        liste=asm_orf.split("_")
        str1 = "_"
        contig=str1.join(liste[:-1])
        if contig in anvio_contig_list:
            if contig in contig2taxo :
                contig2taxo[contig].append(gene_taxo)
            else :
                contig2taxo[contig]=[gene_taxo]

    print(f"length contig2taxo: {len(contig2taxo)}")
    return contig2taxo

def get_contig_taxonomy(list_taxo, ratio):
    """ give a consensus taxonomy from a list of taxonomy  
    :param list_taxo: list of taxonomy
    :param ratio: threshold of taxonomy occurrence
    :return: consensus taxonomy as a LIST
    """
    taxo_contig = []   
    domain_list = []
    phylum_list = []
    class_list = []
    order_list = []
    family_list = []
    genus_list = []
    species_list = []
    
    for gene_taxo in list_taxo :
    
        liste = gene_taxo.split(";")
        domain, phylum, classe, order, family, genus, species = liste[0:7]
        
        domain_list.append(domain)
        phylum_list.append(phylum)
        class_list.append(classe)
        order_list.append(order)
        family_list.append(family)
        genus_list.append(genus)
        species_list.append(species)
        
    domain_counts = Counter(domain_list)
    phylum_counts = Counter(phylum_list)
    class_counts = Counter(class_list)
    order_counts = Counter(order_list)
    family_counts = Counter(family_list)
    genus_counts = Counter(genus_list)
    species_counts = Counter(species_list)
    
    max_domain = max(domain_counts, key=domain_counts.get)
    max_domain_count = domain_counts[max_domain]
    
    max_phylum = max(phylum_counts, key=phylum_counts.get)
    max_phylum_count = phylum_counts[max_phylum]
        
    max_class = max(class_counts, key=class_counts.get)
    max_class_count = class_counts[max_class]
      
    max_order = max(order_counts, key=order_counts.get)
    max_order_count = order_counts[max_order]
        
    max_family = max(family_counts, key=family_counts.get)
    max_family_count = family_counts[max_family]
        
    max_genus = max(genus_counts, key=genus_counts.get)
    max_genus_count = genus_counts[max_genus]
        
    max_species = max(species_counts, key=species_counts.get)
    max_species_count = species_counts[max_species]
        
        
    if float(max_domain_count)/ float(len(list_taxo)) >= ratio:
        contig_domain = max_domain
    else:
        contig_domain = "d__"
        
    if float(max_phylum_count)/ float(len(list_taxo)) >= ratio:
        contig_phylum = max_phylum
    else:
        contig_phylum = "p__"
        
    if float(max_class_count)/ float(len(list_taxo)) >= ratio:
        contig_class = max_class
    else:
        contig_class = "c__"
        
    if float(max_order_count)/ float(len(list_taxo)) >= ratio:
        contig_order = max_order
    else:
        contig_order = "o__"
        
    if float(max_family_count)/ float(len(list_taxo)) >= ratio:
        contig_family = max_family
    else:
        contig_family = "f__"
        
    if float(max_genus_count)/ float(len(list_taxo)) >= ratio:
        contig_genus = max_genus
    else:
        contig_genus = "g__"
        
    if float(max_species_count)/ float(len(list_taxo)) >= ratio:
        contig_species = max_species
    else:
        contig_species = "s__"
        
    taxo_contig=[contig_domain, contig_phylum, contig_class, contig_order, contig_family, contig_genus, contig_species]
    
    return taxo_contig

def get_scaffold_infos(contig2taxo, ratio, contigs_infos, contigs_scgs, contigs_taxo, split_coverage, split_bin):
    """ give all the contig informations in one tuple
    :param contig2taxo: dictionary for contig to gene taxonomy as list
    :param ratio: threshold of taxonomy occurrence
    :param contigs_infos: list of tuples resulting from sql_infos_query statement
    :param contigs_scgs: list of tuples resulting from sql_scgs_query statement
    :param contigs_taxo: list of tuples resulting from sql_taxo_query statement
    :param split_coverage: list of tuples resulting from sql_coverage_query statement
    :param split_bin: list of tuples resulting from sql_bin_query statement
    :return: dictionary for contig to contigs informations
    """

    scaffold_infos = defaultdict(lambda: {
    'gtdb_domain': None, 'gtdb_phylum': None, 'gtdb_class': None, 'gtdb_order': None, 'gtdb_family': None, 'gtdb_genus': None, 'gtdb_species': None,
    'length': None, 'gc_content': None, 'nb_splits': None, 'nb_genes': None, 'coverage': None, 'bin_name': 'Unbinned', 'anvio_genus': None,
    'Archaea_76': [], 'Bacteria_71': [], 'Protista_83': [], 'split_discrepancies': []
    })
    
    # get contig gtdb taxonomy from a list of gene taxonomy
    print('get contigs GTDB taxonomy...')
    gtdb_ranks = ('gtdb_domain', 'gtdb_phylum', 'gtdb_class', 'gtdb_order', 'gtdb_family', 'gtdb_genus', 'gtdb_species')

    for contig in contig2taxo :
        taxo_contig = get_contig_taxonomy(contig2taxo[contig], ratio)

        scaffold_infos[contig].update({gtdb_ranks[i]: taxo_contig[i] for i in range(len(taxo_contig))})

    # get contigs infos from anvio's contig_db
    print('get contigs infos from anvio contig_db...')
    #contigs_infos = get_contigs_info(contigs_db, sql_infos_query)
    vals = ('contig_id', 'length', 'gc_content', 'nb_splits', 'nb_genes')

    for contig in contigs_infos:
        contig_name = contig[0]
        scaffold_infos[contig_name].update({vals[i]: contig[i] for i in range(1, len(contig))})

    # get anvio hmms Single Copy Genes to compute completeness
    print('get anvio hmms Single Copy Genes...')
    #contigs_scgs = get_contigs_info(contigs_db, sql_scgs_query)
    for contig in contigs_scgs: 
        contig_name = contig[0].split("_split")[0]
        source = contig[1]
        gene_caller_id = contig[2]
        gene_name = contig[3]
        gene_id = gene_name + " ("+ str(gene_caller_id) +")"
    
        scaffold_infos[contig_name].setdefault(source, []).append(gene_id)

    # get anvio taxonomy (genus rank)
    print('get anvio taxonomy...')
    anvio_contig_taxo = defaultdict(list)
    #contigs_taxo = get_contigs_info(contigs_db, sql_taxo_query)
    for contig in contigs_taxo:
        contig_name = contig[0]
        genus = contig[2]
    
        if contig_name not in anvio_contig_taxo:
            anvio_contig_taxo[contig_name] = [genus]
        else:
            anvio_contig_taxo[contig_name].append(genus)
    
    for contig, taxo in anvio_contig_taxo.items():
        genus_counts = Counter(taxo)        
        majority_genus = max(genus_counts, key=genus_counts.get)

        #scaffold_infos[contig]['anvio_genus'] = majority_genus
        scaffold_infos[contig].update({'anvio_genus': majority_genus})
            
    # get coverage and Bin_id from one collection in anvio PROFILE_db
    print(f'get coverage and Bin_id from anvio collection {collection}...')
    split_infos = defaultdict(lambda: {'coverage': None,'bin_name': 'Unbinned'})
    contig_splits = defaultdict(list)
    split_outliers = {}

    #split_coverage = get_contigs_info(profile_db, sql_coverage_query)
    for contig in split_coverage:
        split_name = contig[0]
        coverage = float(contig[1])
        split_infos[split_name].update({'coverage': coverage})

    #split_bin = get_bins_info(profile_db, sql_bin_query, collection)
    for contig in split_bin:
        split_name = contig[0]
        bin_name = contig[1]
        split_infos[split_name].update({'bin_name': bin_name})

    for split_name, infos in split_infos.items():
        contig_name = split_name.split("_split")[0]
        coverage = infos.get('coverage', 0.0)
        bin_name = infos.get('bin_name', 'Unbinned')
    
        contig_splits[contig_name].append((split_name, coverage, bin_name))

    for contig, splits in contig_splits.items():
        coverage_counts = Counter([split[1] for split in splits])
        bin_counts = Counter([split[2] for split in splits])
        
        majority_coverage = max(coverage_counts, key=coverage_counts.get)
        majority_bin = max(bin_counts, key=bin_counts.get)
        
        scaffold_infos[contig].update({'coverage': majority_coverage, 'bin_name': majority_bin})             
                
        logging.basicConfig(filename="warnings.log", level=logging.WARNING, format="%(levelname)s: %(message)s")

        for split_name, coverage, bin_name in splits:
            outlier_type = []
            if coverage != majority_coverage:
                outlier_type.append("COV")
                
            if bin_name != majority_bin:
                outlier_type.append("BIN")
                
            if outlier_type:
               # logging.warning(f"{split_name} is an outlier for {outlier_type[0]}")
                logging.warning("{} is an outlier for {}".format(split_name, outlier_type[0]))

                split_outliers[split_name] = {
                    'coverage': coverage,
                    'bin_name': bin_name,
                    'outliers': ",".join(outlier_type)  
                }
                
                split_nb = "split"+split_name.split("_split")[1]
                
                scaffold_infos[contig].setdefault('split_discrepancies', []).append(split_nb + " " + ",".join(outlier_type))

    print(f"length scaffold_infos: {len(scaffold_infos)}")
    return split_outliers,scaffold_infos

def running_anvi_summarize(collection, profileDb_filename, contigDb_filename):
    """ run anvi-summarize  
    :param collecion: anvio collection to summarize 
    :param profileDb_filename: anvio PROFILE database
    :param contigDb_filename: anvio contigs database
    :return: 
    """
    anvio_summary_dir = wdir + "/MAGnet/ANVIO-SUMMARY"
    cmd = 'source activate '+anvioVersion+' && anvi-summarize --report-aa-seqs-for-gene-calls -p '+profileDb_filename+' -c '+contigDb_filename+' -o '+anvio_summary_dir+' -C '+collection+' >/dev/null 2>&1'
    print(cmd)
    status = os.system(cmd)
    #print('status: '+str(status))

    if not status == 0:
        sys.exit("something went wrong with anvi-summarize, exit")

def get_project_sample_names(contigs_filename) :
    """ extract  Project and Sample name from contigs_names
    :param contigs_filename: fasta file contanning contigs renamed by anvio
    :return: project and sample names in text format
    """
    file = open(contigs_filename,'r')
    for line in file :
        if re.match('>',line) :
            name = line.rstrip().replace('>','')
            project,sample,contig = name.split('__')
            break
        else:
            continue
    file.close()
    return project,sample

def suggested_name_function(anvio_lineage, gtdb_lineage, project, sample, binName2count):
    """ give a uniq name to each bin
    :param anvio_lineage: bin anvio taxonomy
    :param gtdb_lineage: bin gtdb taxonomy 
    :param project: project name
    :param sample: sample name
    :return: bin name in text format
    """
    name = 'Unknown'

    if gtdb_lineage and gtdb_lineage != 'NULL':
        liste = gtdb_lineage.split(';')
        if re.match(r's\_\_',liste[-1]) :
            if liste[-1] == 's__' : # if gtdb species name is empty, look for genus, family, order, class, phylum, domain name...
                for taxa in reversed(liste) :
                    if taxa.split('__')[1] != '' :
                        name = taxa.split('__')[1]
                        break
                    else:
                        continue
            else: # if has a gtdb species name then concatenate genus and species name
                name = liste[-1].split('__')[1].replace(' ','_')
        else: # no species s__ field...
            name = liste[-1].split('__')[1]
        name = name.replace('-','_')
    else:
        liste = anvio_lineage.split(';')
        for taxa in reversed(liste) :
            if taxa.strip() != '':
                name = taxa.strip()
                break
            else:
                continue
        name = name.replace(' ','_')
        name = name.replace('-','_')

    binName = name+'__'+project+'__'+sample
    # check name redundancy
    if binName in binName2count :
        binName2count[binName] += 1
        nb = str(binName2count[binName])
        binName = name+'_'+nb+'__'+project+'__'+sample
    else:
        binName2count[binName] = 1

    return binName

def link2bins_contigs_files(anvi_summarize_directory, binned_list):
    bins_dir = wdir+ '/MAGnet/bins'
    if os.path.exists(bins_dir) :
        sys.exit(bins_dir+' already exist, please remove it first')
    os.mkdir(bins_dir)
    
    for root, dirs, files in os.walk(anvi_summarize_directory+ '/bin_by_bin', topdown = False):
        for binName in dirs:
            fasta_filename = anvi_summarize_directory+'/bin_by_bin/'+ binName +'/'+ binName + '-contigs.fa'
            os.symlink(fasta_filename,bins_dir +'/'+ binName +'.fna')

    #creating an unbinned file
    unbinned_list = list()
    for record in SeqIO.parse(anvio_contig_filename,'fasta'):
        if record.id not in binned_list:
            unbinned_list.append(record)
            
            unbinned_filename = bins_dir+ '/Unbinned.fna'
            SeqIO.write(unbinned_list,unbinned_filename,'fasta')

#### choisir la version de GTDB-tk
def runningGTDBtk(gtdbtk_dir,bins_dir,cpu) :
    """ run GTDB-Tk on each bin
    :param gtdbtk_dir: directory to store bins and results
    :param bins_dir: directory that contains links to the fasta bin files 
    :param cpu: number of threads 
    :return:
    """
    os.mkdir(gtdbtk_dir)
    os.mkdir(gtdbtk_dir+'/'+'bins')
    os.mkdir(gtdbtk_dir+'/'+'output')
    
    for root, dirs, files in os.walk(bins_dir, topdown = False):
        for filename in files :
            if re.search(r'.fna',filename) :
                print(filename)
                binName = filename.replace('.fna','')
                if binName == 'Unbinned' :
                    continue
                print(root +'/'+ filename)
                os.symlink(root +'/'+ filename, gtdbtk_dir +'/bins/'+ binName +'.fna')
    #cmd = 'source activate gtdbtk-2.4.0 && gtdbtk classify_wf  --skip_ani_screen --cpus ' +cpu+' --genome_dir ' +gtdbtk_dir+ '/bins' + ' --out_dir ' +gtdbtk_dir+ '/output > ' +gtdbtk_dir+ '/gtdbtk.log 2>&1'
    cmd = 'source activate gtdbtk-1.7.0 && gtdbtk classify_wf  --cpus '+cpu+' --genome_dir '+gtdbtk_dir+'/'+'bins'+' --out_dir '+gtdbtk_dir+'/'+'output > '+gtdbtk_dir+'/'+'gtdbtk.log 2>&1'
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status)+'\n')
    if not status == 0 :
        sys.exit('something went wrong with gtdbtk classify_wf, exit.')

def runningCheckM2(checkm_dir,bins_dir,cpu) :
    """ run CheckM2 on each bin
    :param checkm_dir: directory to store bins and results
    :param bins_dir: directory that contains links to the fasta bin files 
    :param cpu: number of threads 
    :return:
    """
    os.mkdir(checkm_dir)
    os.mkdir(checkm_dir+'/'+'bins')
    os.mkdir(checkm_dir+'/'+'output')
    
    for root, dirs, files in os.walk(bins_dir, topdown = False):
        for filename in files :
            if re.search(r'.fna',filename) :
                print(filename)
                binName = filename.replace('.fna','')
                if binName == 'Unbinned' :
                    continue
                print(root+'/'+filename)
                os.symlink(root+'/'+filename,checkm_dir+'/'+'bins'+'/'+binName+'.fna')
    cmd = 'source activate checkm2 && checkm2 predict --threads ' +cpu+ ' --input ' +checkm_dir+ '/bins' + ' --output-directory ' +checkm_dir+ '/output > ' +checkm_dir+ '/checkm2.log 2>&1'
    print(cmd)
    status = os.system(cmd)
    print('status: '+str(status)+'\n')
    if not status == 0 :
        sys.exit('something went wrong with checkm lineage_wf, exit.')

def get_bins_infos(anvio_summary_dir, gtdbtk_output, checkm2_output) :
    """ give all the bin informations in one tuple
    :param anvio_summary_dir: directory containing results of anvi-summarize
    :param gtdb-tk_output: directory containing results of GTDB-Tk
    :param checkm2_output: directory containing results of CheckM2
    :return: dictionary for bin to all bins informations
    """
    bins_infos = defaultdict(lambda: {
    'name': None, 'coverage': None, 'length': None, 'contig_nb': None, 'N50': None, 'gc': None, 'completeness': None, 'contamination': None, 'taxonomy': None,
    'gtdb_taxonomy': None, 'gtdb_fastani_reference': None, 'gtdb_classification_method': None,
    'checkm2_completeness_model': None, 'checkm2_translation_table': None, 'checkm2_coding_density': None, 'checkm2_completeness': None, 'checkm2_contamination': None, 'anvio_bin': '1'
})

    with open(anvio_summary_dir+'/bins_across_samples/mean_coverage.txt') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        for line in reader:
            bin_name = line[0]
            coverage = line[1]
        
            bins_infos[bin_name].update({'coverage': coverage})

    with open(anvio_summary_dir+'/bins_summary.txt') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        liste_taxo = []
        for line in reader:
            bin_name = line[0]
            length = int(line[1])
            contig_nb = int(line[2])
            N50 = int(line[3])
            gc = float(line[4])/100
            completeness = float(line[5])/100
            contamination = float(line[6])/100
            
            liste_taxo = [line[i] for i in range(7, 14) if i < len(line) and line[i].strip()]
            taxonomy = ";".join(liste_taxo)
           
            bins_infos[bin_name].update({
                'bins_name': bin_name,
                'length': length,
                'contig_nb': contig_nb,
                'N50': N50,
                'gc': gc,
                'completeness': completeness,
                'contamination': contamination,
                'taxonomy': taxonomy
            })
    
    with open(gtdbtk_output+'/classify/gtdbtk.bac120.summary.tsv') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        for line in reader:
            bin_name = line[0]
            taxonomy = line[1]
            reference = line[2]
            method = line[13]
    
            bins_infos[bin_name].update({
                'gtdb_taxonomy': taxonomy,
                'gtdb_fastani_reference': reference,
                'gtdb_classification_method': method
            })

    with open(checkm2_output+'/quality_report.tsv') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        vals = ('bin_name', 'checkm2_completeness', 'checkm2_contamination', 'checkm2_completeness_model', 'checkm2_translation_table', 'checkm2_coding_density')
        for line in reader:
            bin_name = line[0]
            bins_infos[bin_name].update({vals[i]: line[i] for i in range(1, len(vals))})
         
    return bins_infos       




if __name__ == "__main__":
    #/env/cns/bigtmp2/stephanie/test_MAGnet_120225/gMMP_SC1_AE1

    parser = argparse.ArgumentParser(description="Run the MAGnet pipeline")    
    parser.add_argument("-config", help="the path to the info.json file created by the metagenomic pipeline")
    parser.add_argument("-collection", help="the name of the anvio collection that contain bins to clean")
    parser.add_argument("-cpu",type=int, default=2, help="number of CPUs to use")
    args = parser.parse_args()

    # checking arguments

    if args.config == None :
        sys.exit('the config parameter is empty, please provide one json file')
    else:
        config_file = args.config

    if args.collection == None :
        sys.exit('the collection parameter is empty, please provide one collection')
    else:
        collection = args.collection

    if args.cpu < 2 :
        cpu = str(2)
    else:
        cpu = str(args.cpu)

    with open(config_file) as f:
        data = json.load(f)

    wdir = data['directory']
    profileDb_filename = data['anvio_profileDb_filename']
    contigsDb_filename = data['anvio_contigDb_filename']
    assembly_protein_filename = data['assembly_protein_filename']
    anvio_contig_filename = data['anvio_contig_filename']
    anvioVersion = data['anvio_version']

    print('\n')
    print('##############')
    print('# parameters #')
    print('##############')
    print('\n')
    print('anvio PROFILE.db: '+profileDb_filename)
    print('anvio contigs.db: '+contigsDb_filename)
    print('Collection name: '+collection)
    print('Working directory: '+wdir)
    print('Number of CPUs: '+str(cpu))

    print('\n')
    print('############')
    print('# pipeline #')
    print('############')
    print()

    MAGnet_directory = wdir+ "/MAGnet"
    if os.path.exists(MAGnet_directory) :
        sys.exit(MAGnet_directory+" already exist, please rename it first")
    os.mkdir(MAGnet_directory)

    print('\n')
    print('##############################')
    print('# Step 1:get scaffolds infos #')
    print('##############################')

    print('get GTDB taxonomy...')
    gtdb_path = "/env/cns/proj/agc/bank/GTDB/release220/"

    print(f"path to the GTDB release: {gtdb_path}")
    archaea_metadata_filename = gtdb_path + "ar53_metadata_r220.tsv"
    bacteria_metadata_filename = gtdb_path + "bac120_metadata_r220.tsv"
    
    access2taxo = merge_dicts(get_GTDB_taxonomy(archaea_metadata_filename), get_GTDB_taxonomy(bacteria_metadata_filename))
    print(f"length access2taxo: {len(access2taxo)}")

    gtdb_protein_filename = gtdb_path + "gtdb_r220_proteins_db.faa"
    #orf2access = get_GTDB_proteins(gtdb_protein_filename)
    
    diamond_path = MAGnet_directory + "/diamond"
    if os.path.exists(diamond_path) :
        shutil.rmtree(diamond_path)
    os.mkdir(diamond_path)

    k=3
    
    print(diamond_path)

    start_time = time.time()
    print(f"run diamond blastp...")
    gtdb_protein_db = gtdb_path + "gtdb_r220_proteins_db.faa"
    #run_diamond_blastp(diamond_path, gtdb_protein_db, assembly_protein_filename, cpu)
    end_time = time.time()
    print(f"diamond duration = {end_time - start_time:.2f} seconds")

    asm_orf2best_hits = get_diamond_besthits(diamond_path)

    with open(MAGnet_directory+'/diamond_best_hits.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['asm_orf','bitscore', 'accession', 'GTDB_orf','diamond_besthit'])
        
        for asm_orf in asm_orf2best_hits:
            cpt = 0
            besthit_list=[]
            #print(asm_orf)
            for bitscore in sorted(asm_orf2best_hits[asm_orf],reverse=True):

                #print('\t'+str(bitscore))
                for access_orf in asm_orf2best_hits[asm_orf][bitscore]:
                    cpt += 1
                    access,orf = access_orf.split("__") 
                    taxo = access2taxo[access]
                    
                if cpt <= 3 :
                    besthit_list.append(taxo)
                    writer.writerow([asm_orf+'\t'+str(bitscore)+'\t'+access+'\t'+orf+'\t'+taxo])
                    writer.writerow([asm_orf+'\t'+'consensus'+'\t'+k_best_hits(asm_orf2best_hits[asm_orf],k,access2taxo)])
                    
                else:
                    break
         
 


    anvio_contigs_list=[]
    for record in SeqIO.parse(anvio_contig_filename,'fasta'):
        anvio_contigs_list.append(record.id)
    print(f"number of contigs in anvio: {len(anvio_contigs_list)}")

    contig2taxo = get_gene_taxonomy(anvio_contigs_list, asm_orf2best_hits, k, access2taxo)
    

    #sys.exit()

    print('\n')
    print('############')
    print('get contigs infos from ANVIO databases...') 
    print('\n') 
  
    
    sql_infos_query = """SELECT contigs_basic_info.contig,
        length,
        gc_content,
        num_splits,
        COUNT(genes_in_contigs.contig) 
        FROM contigs_basic_info 
        INNER JOIN genes_in_contigs 
        ON contigs_basic_info.contig = genes_in_contigs.contig 
        GROUP BY genes_in_contigs.contig 
    ;"""

    sql_scgs_query = """SELECT split,
        source,
        hmm_hits.gene_callers_id,
        gene_name,
        gene_hmm_id
        FROM genes_in_splits
        INNER JOIN hmm_hits
        ON genes_in_splits.gene_callers_id = hmm_hits.gene_callers_id
    ;"""
    
    sql_taxo_query = """SELECT contig,
        genes_taxonomy.gene_callers_id,
        taxon_names.t_genus
        FROM genes_taxonomy
        INNER JOIN genes_in_contigs
        ON genes_taxonomy.gene_callers_id = genes_in_contigs.gene_callers_id
        INNER JOIN taxon_names
        ON genes_taxonomy.taxon_id = taxon_names.taxon_id
    ;"""
    
    sql_coverage_query = """SELECT item,
        value
        FROM mean_coverage_contigs
    ;"""

    sql_bin_query = """SELECT split,
        bin_name
        FROM collections_of_splits
        WHERE collection_name = ? 
    ;"""

    contigs_infos = get_contigs_info(contigsDb_filename, sql_infos_query)
    contigs_scgs = get_contigs_info(contigsDb_filename, sql_scgs_query)
    contigs_taxo = get_contigs_info(contigsDb_filename, sql_taxo_query)
    split_coverage = get_contigs_info(profileDb_filename, sql_coverage_query)
    split_bin = get_bins_info(profileDb_filename, sql_bin_query, collection)
    
    ratio=0.5
    split_outliers,scaffold_infos = get_scaffold_infos(contig2taxo, ratio, contigs_infos, contigs_scgs, contigs_taxo, split_coverage, split_bin)
    
    print('\n')
    print('##############################')
    print('# Step 2:get bins infos #')
    print('##############################')

    print('\n')
    print('running_anvi_summarize...')
    print('\n')

    running_anvi_summarize(collection, profileDb_filename, contigsDb_filename)

    print('create bins directory...')
    anvio_summary_dir = wdir + '/MAGnet/ANVIO-SUMMARY'

    binned_list = []
    for scaffold, infos in scaffold_infos.items():
        if infos.get('bin_name') != 'Unbinned':
            binned_list.append(scaffold)
    
    link2bins_contigs_files(anvio_summary_dir, binned_list)
    
    print('\n')
    print('running GTDBtk on bins...')
    print('\n')
    start_time = time.time()
    bins_dir = wdir+ '/MAGnet/bins'
    gtdbtk_dir = wdir+ '/MAGnet/GTDB-tk'
    runningGTDBtk(gtdbtk_dir,bins_dir,cpu)
    end_time = time.time()
    print(f"GTDBtk duration = {end_time - start_time:.2f} seconds")

    print('runningCheckM2 on bins...')
    print('\n')
    start_time = time.time()
    checkm_dir = wdir+ '/MAGnet/CheckM2'
    runningCheckM2(checkm_dir,bins_dir,cpu)
    end_time = time.time()
    print(f"CheckM2 duration = {end_time - start_time:.2f} seconds")

    print('############')
    print('get bins infos from ANVIO, GTDB-tk and CheckM2...')
    print('\n')
    anvio_summary_dir = wdir + '/MAGnet/ANVIO-SUMMARY'
    gtdbtk_output = gtdbtk_dir+ '/output'
    checkm2_output = checkm_dir+ '/output'
    bins_infos = get_bins_infos(anvio_summary_dir, gtdbtk_output, checkm2_output)
     
    project,sample = get_project_sample_names(anvio_contig_filename)
    
    binName2count={}
    for bin_name, bin_data in bins_infos.items():
        anvio_lineage = str(bin_data['taxonomy'])
        gtdb_lineage = str(bin_data['gtdb_taxonomy']) if bin_data["gtdb_taxonomy"] else "NULL"
         
        binName = suggested_name_function(anvio_lineage, gtdb_lineage, project, sample, binName2count)
    
        bins_infos[bin_name].update({'name': binName})

    print('\n')
    print('##############################')
    print('# Step 3:create MAGnet.db    #')
    print('##############################')
    print('\n')

    # Creation of the MAGnet database and the scaffolds table
    print('create MAGnet database and the scaffolds table...')
    magnet_db_file = wdir + '/MAGnet' + '/MAGNET.db'
    
    create_scaffolds_statements = [""" CREATE TABLE IF NOT EXISTS scaffolds (
                scaffold_name PRIMARY KEY,
                gtdb_domain text,
                gtdb_phylum text,
                gtdb_class text,
                gtdb_order text,
                gtdb_family text,
                gtdb_genus text,
                gtdb_species text,
                anvio_length integer,
                anvio_gc real,
                anvio_nb_splits integer,
                anvio_nb_genes integer,
                anvio_coverage real,
                anvio_updated_bin_id text NOT NULL,
                anvio_bin_id text NOT NULL,
                anvio_genus text,
                anvio_scgs_Archaea_76 text,
                anvio_scgs_Bacteria_71 text,
                anvio_scgs_Protista_83 text,
                split_discrepancies text,
                FOREIGN KEY (anvio_bin_id) REFERENCES bins (anvio_id)
        );"""]

    scaffolds_keys_order = [
        'gtdb_domain', 'gtdb_phylum', 'gtdb_class', 'gtdb_order', 'gtdb_family', 
        'gtdb_genus', 'gtdb_species', 'length', 'gc_content', 'nb_splits', 'nb_genes', 'coverage', 
        'bin_name', 'bin_name', 'anvio_genus', 'Archaea_76', 'Bacteria_71', 'Protista_83', 'split_discrepancies'
    ]

    sql_add_scaffold = ''' INSERT INTO scaffolds(scaffold_name,gtdb_domain,gtdb_phylum,gtdb_class,gtdb_order,gtdb_family,gtdb_genus,gtdb_species,
        anvio_length,anvio_gc,anvio_nb_splits,anvio_nb_genes,anvio_coverage,anvio_updated_bin_id,anvio_bin_id,anvio_genus,
        anvio_scgs_Archaea_76,anvio_scgs_Bacteria_71,anvio_scgs_Protista_83,split_discrepancies)
            VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'''

    create_table(magnet_db_file, create_scaffolds_statements)
    print('fill scaffolds table...')
    list_scaffolds = data2insert(scaffold_infos, scaffolds_keys_order)
    populate_table(magnet_db_file, sql_add_scaffold, list_scaffolds)

    print('create bins table...')
    create_bins_statements = [""" CREATE TABLE IF NOT EXISTS bins (
                anvio_id text PRIMARY KEY,
                name text,
                anvio_coverage real,
                anvio_length integer,
                anvio_contig_nb integer,
                anvio_N50 integer,
                anvio_gc real,
                anvio_completeness real,
                anvio_contamination real,
                anvio_taxonomy text,
                gtdb_taxonomy text,
                gtdb_fastani_reference text,
                gtdb_classification_method text,
                checkm2_completeness_model text,
                checkm2_translation_table integer,
                checkm2_coding_density real,
                checkm2_completeness real,
                checkm2_contamination real,
                anvio_bin integer NOT NULL
        );"""]

    bins_keys_order = [
        'name', 'coverage', 'length', 'contig_nb', 'N50', 'gc', 'completeness', 'contamination', 'taxonomy',
        'gtdb_taxonomy', 'gtdb_fastani_reference', 'gtdb_classification_method',
        'checkm2_completeness_model', 'checkm2_translation_table', 'checkm2_coding_density', 'checkm2_completeness', 'checkm2_contamination', 'anvio_bin'
    ]

    sql_add_bin = """ INSERT INTO bins (anvio_id,name,anvio_coverage,anvio_length,anvio_contig_nb,anvio_N50,anvio_gc,
            anvio_completeness,anvio_contamination,anvio_taxonomy,gtdb_taxonomy,gtdb_fastani_reference,gtdb_classification_method,
            checkm2_completeness_model,checkm2_translation_table,checkm2_coding_density,checkm2_completeness,checkm2_contamination,anvio_bin)
                  VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"""
    
    print('\n')
    create_table(magnet_db_file, create_bins_statements)
    print('fill bins table...')
    list_bins = data2insert(bins_infos, bins_keys_order)
    populate_table(magnet_db_file, sql_add_bin, list_bins)

