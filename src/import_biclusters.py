#!/usr/bin/env python

"""
This import script takes Matt Wall's JSON scripts and stuffs them into
a MySQL database.
"""
import os
import sys
import json
import mysql.connector
import traceback

DATABASE_NAME = 'mm_api'
DATABASE_USER = 'egrin2'
DATABASE_PASS = 'egrin2'

# quick mapping of roles to database id
MUTATION_TF_ROLES = { 'down_regulates': 1, 'up_regulates': 2 }
TF_BC_ROLES = { 'activates': 1, 'represses': 2 }

def collect_genes(data, genes):
    for bicluster_id, bcdata in data.items():
        genes.update(bcdata['genes'])

def collect_tfs(data, tfs):
    for bicluster_id, bcdata in data.items():
        edge1 = bcdata['edge_1']
        edge2 = bcdata['edge_2']
        tfs.add(edge1[2])
        tfs.add(edge2[0])

def collect_mutations(data, mutations):
    for bicluster_id, bcdata in data.items():
        edge1 = bcdata['edge_1']
        mutations.add(edge1[0])


def populate_database(directory, conn):
    genes = set()
    tfs = set()
    mutations = set()
    mutation_tf_roles = set()
    tf_bc_roles = set()
    files = [f for f in os.listdir(directory) if f.endswith('json')]
    for f in files:
        path = os.path.join(sys.argv[1], f)
        with open(path) as infile:
            data = json.load(infile)
            collect_genes(data, genes)
            collect_tfs(data, tfs)
            collect_mutations(data, mutations)
    cursor = conn.cursor()
    try:
        for gene in sorted(genes):
            cursor.execute('insert into genes (name) values (%s)', [gene])
        for tf in sorted(tfs):
            cursor.execute('insert into tfs (name) values (%s)', [tf])
        for mutation in sorted(mutations):
            cursor.execute('insert into mutations (name) values (%s)', [mutation])
    finally:
        cursor.close()


def get_maps(conn):
    gene_map = {}
    mutation_map = {}
    tf_map = {}
    cursor = conn.cursor()
    try:
        cursor.execute('select id,name from genes')
        gene_map = {name: pk for pk, name in cursor.fetchall()}
        cursor.execute('select id,name from mutations')
        mutation_map = {name: pk for pk, name in cursor.fetchall()}
        cursor.execute('select id,name from tfs')
        tf_map = {name: pk for pk, name in cursor.fetchall()}
    finally:
        cursor.close()
    return gene_map, mutation_map, tf_map


def insert_biclusters(directory, conn):
    gene_map, mutation_map, tf_map = get_maps(conn)
    files = [f for f in os.listdir(directory) if f.endswith('json')]
    cursor = conn.cursor()
    try:
        # insert biclusters to establish them in the database
        for f in files:
            path = os.path.join(sys.argv[1], f)
            with open(path) as infile:
                data = json.load(infile)
                for bicluster_id, bcdata in data.items():
                    mutation, bcnum = bicluster_id.split('_')
                    mutation_id = mutation_map[mutation]
                    cursor.execute('insert into biclusters (mutation_id,name) values (%s,%s)',
                                   [mutation_id, bicluster_id])
                    bc_pk = cursor.lastrowid
                    for g in bcdata['genes']:
                        cursor.execute('insert into bicluster_genes (bicluster_id,gene_id) values (%s,%s)',
                                       [bc_pk, gene_map[g]])
                    edge1 = bcdata['edge_1']
                    cursor.execute('insert into bc_mutation_tf (bicluster_id,mutation_id,tf_id,role) values (%s,%s,%s,%s)',
                                   [bc_pk, mutation_map[edge1[0]], tf_map[edge1[2]],
                                    MUTATION_TF_ROLES[edge1[1]]])

                    edge2 = bcdata['edge_2']
                    cursor.execute('insert into bc_tf (bicluster_id,tf_id,role) values (%s,%s,%s)',
                                   [bc_pk, tf_map[edge2[0]], TF_BC_ROLES[edge2[1]]])
        conn.commit()
    except:
        traceback.print_exc()
    finally:
        cursor.close()

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print("usage: import_biclusters.py <directory>")
    else:
        conn = mysql.connector.connect(user=DATABASE_USER,
                                       password=DATABASE_PASS,
                                       database=DATABASE_NAME)
        """
        try:
            populate_database(sys.argv[1], conn)
            conn.commit()
        except:
            traceback.print_exc()
            conn.rollback()
        finally:
            conn.close()
        """
        insert_biclusters(sys.argv[1], conn)
