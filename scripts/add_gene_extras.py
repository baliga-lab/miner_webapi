#!/usr/bin/env python3

import os
import requests
import xml.etree.ElementTree as ET
import re
import MySQLdb
from collections import defaultdict

DBNAME = 'gbm_api'


def dbconn():
    return MySQLdb.connect(host="localhost", user="admin", passwd="password",
                           db=DBNAME)


def get_gene_extra(ensembl_id, query_description=True, query_function=True):
    """Retrieve extra information from EnsEMBL. Takes a long time, don't call this too often"""
    if query_description:
        r = requests.get('https://rest.ensembl.org/lookup/id/%s?content-type=application/json;expand=1' % ensembl_id)
        if r.status_code == 200:
            ensdata = r.json()
            try:
                description = ensdata['description']
            except KeyError:
                description = None
        else:
            description = None
    else:
        description = None

    r = requests.get('https://rest.ensembl.org/xrefs/id/' + ensembl_id +
                     '?content-type=application/json;external_db=uniprot%;all_levels=1')
    if r.status_code == 200:
        xrefdata = r.json()
        uniprot_ids = [entry['primary_id'] for entry in xrefdata]
        if len(uniprot_ids) > 0:
            uniprot_id = uniprot_ids[0]
        else:
            uniprot_id = None
    else:
        uniprot_id = None

    if query_function and uniprot_id is not None:
        # retrieve the function information from UniProt
        r = requests.get('https://www.uniprot.org/uniprot/%s.xml' % uniprot_id)
        doc = ET.fromstring(r.text)
        function = '-'
        for child in doc:
            localname = re.sub(r'{.*}', '', child.tag)
            if localname == 'entry':
                for c in child:
                    localname = re.sub(r'{.*}', '', c.tag)
                    if localname == 'comment' and c.attrib['type'] == 'function':
                        function = ""
                        for node in c:
                            function += node.text
    else:
        function = None

    return uniprot_id, description, function


if __name__ == '__main__':
    conn = dbconn()
    cur = conn.cursor()
    try:
        cur.execute('select ensembl_id from genes where uniprot_id is null and ensembl_id is not null')
        ensembl_ids = [row[0] for row in cur.fetchall() if row[0] is not None]
        for index, ensembl_id in enumerate(ensembl_ids):
            print("Update information for '%s' (%d of %d)" % (ensembl_id, index + 1, len(ensembl_ids)))
            uniprot_id, description, function = get_gene_extra(ensembl_id)
            cur.execute('update genes set uniprot_id=%s,uniprot_function=%s,ens_description=%s where ensembl_id=%s',
            [uniprot_id, description, function, ensembl_id])
            conn.commit()
    finally:
        cur.close()
