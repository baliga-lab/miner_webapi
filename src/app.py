#!/usr/bin/env python

import logging
import json
import os
import traceback
from collections import defaultdict
import mysql.connector

# for mockup
import random

from flask import Flask, Response, url_for, redirect, render_template, request, session, flash, jsonify
import flask
from sqlalchemy import and_
# because we have an API, we need to allow cross-origin here
from flask_cors import CORS, cross_origin

import requests
import xml.etree.ElementTree as ET
import re
import datasets

app = Flask(__name__)
CORS(app)

app.config.from_envvar('APP_SETTINGS')

MUTATION_REGULATOR_ROLES = { 1: 'down-regulates', 2: 'up-regulates'}
REGULATOR_REGULON_ROLES = { 1: 'activates', 2: 'represses' }


def dbconn():
    return mysql.connector.connect(user=app.config['DATABASE_USER'],
                                   password=app.config['DATABASE_PASSWORD'],
                                   database=app.config['DATABASE_NAME'])

@app.route('/regulon/<regulon>')
def regulon_info(regulon):
    """all genes in the specified bicluster"""
    conn = dbconn()
    cursor = conn.cursor(buffered=True)
    try:
        cursor.execute("select id,cox_hazard_ratio from regulons where name=%s", [regulon])
        print("REGULON: ", regulon)
        pk, hazard_ratio = cursor.fetchone()
        cursor.execute("select count(*) from cm_flows where regulon_id=%s", [pk])
        num_causalflows = cursor.fetchone()[0]

        # transcription factor -> bicluster
        cursor.execute("""select tfs.ensembl_id,tfs.preferred,role_id
from regulons bc join regulon_regulator bt on bc.id=bt.regulon_id
join genes as tfs on tfs.id=bt.regulator_id
where bc.name=%s""", [regulon])
        regulon_regulators = [
            {
                "regulator": tf,
                "regulator_preferred": tf_preferred if tf_preferred is not None else tf,
                "role": REGULATOR_REGULON_ROLES[role]
            }
            for tf, tf_preferred, role in cursor.fetchall()]

        # regulon genes
        cursor.execute("""select g.preferred
from regulons bc join regulon_genes bcg on bc.id=bcg.regulon_id
join genes g on g.id=bcg.gene_id where bc.name=%s""",
                       [regulon])
        genes = [row[0] for row in cursor.fetchall()]

        # regulon drugs
        cursor.execute("""select d.name
from regulons bc join drug_regulons rd on rd.regulon_id=bc.id
join drugs d on rd.drug_id=d.id where bc.name=%s""",
                       [regulon])
        drugs = [row[0] for row in cursor.fetchall()]

        # transcriptional programs
        cursor.execute("""select tp.name
from regulons reg join regulon_programs rp on reg.id=rp.regulon_id
join trans_programs tp on tp.id=rp.program_id where reg.name=%s""", [regulon])
        programs = [row[0] for row in cursor.fetchall()]

        return jsonify(regulon=regulon,
                       hazard_ratio=hazard_ratio,
                       regulon_regulators=regulon_regulators,
                       genes=genes,
                       drugs=drugs,
                       program=programs[0],  # for now, just return 1 element, it's actually 1:1
                       num_causal_flows=num_causalflows)
    except:
        traceback.print_exc()
    finally:
        cursor.close()
        conn.close()


def _make_causalflow_results(cursor):
    cm_flows = [
        {
            "cmf_id": cmf_id,
            "regulon": regulon,
            "cmf_type": cmf_type,
            "pathway": pathway,
            "mutation_gene_ensembl": mutgene_ensembl,
            "mutation_gene_symbol": mutgene_symbol,
            "mutation": mutation,
            "regulator": regulator_ensembl,
            "regulator_preferred": regulator_symbol if regulator_symbol is not None else regulator_ensembl,
            "regulator_role": REGULATOR_REGULON_ROLES[mutation_regulator_role],
            "mutation_role": MUTATION_REGULATOR_ROLES[1 if regulon_regulator_spearman_r < 0 else 2],
            "regulator_pvalue": regulon_mutation_regulator_pvalue,
            "regulator_spearman_r": regulon_regulator_spearman_r,
            "regulator_spearman_pvalue": regulon_regulator_spearman_pvalue,
            "regulon_t_statistic": regulon_t_statistic,
            "regulon_log10_p_stratification": regulon_log10_p_stratification,
            "fraction_edges_correctly_aligned": fraction_edges_correctly_aligned,
            "fraction_aligned_diffexp_edges": fraction_aligned_diffexp_edges,
            "num_downstream_regulons": num_downstream_regulons,
            "num_diffexp_regulons": num_diffexp_regulons
        }
        for cmf_id, regulon, mutation, cmf_type, pathway,
        mutgene_ensembl, mutgene_symbol,
        regulator_ensembl, regulator_symbol, mutation_regulator_role,
        regulon_mutation_regulator_pvalue,
        regulon_regulator_spearman_r, regulon_regulator_spearman_pvalue,
        regulon_t_statistic, regulon_log10_p_stratification,
        fraction_edges_correctly_aligned, fraction_aligned_diffexp_edges,
        num_downstream_regulons, num_diffexp_regulons
        in cursor.fetchall()]

    return jsonify(cm_flows=cm_flows)


@app.route('/causalflows_for_regulon/<regulon>')
def causalflows_for_regulon(regulon):
    """all causal flows for the specified regulon"""
    conn = dbconn()
    cursor = conn.cursor(buffered=True)
    try:
        cursor.execute("""select cmf_id,r.name,m.name as mutation,cmft.name as cmftype,pw.name as pathway,mg.ensembl_id as mut_ensembl_id,mg.preferred as mut_preferred,tf.ensembl_id as tf_ensembl_id,tf.preferred as tf_preferred,regulon_mutation_regulator_role_id,regulon_mutation_regulator_pvalue,regulon_regulator_spearman_r,regulon_regulator_spearman_pvalue,regulon_t_statistic,regulon_log10_p_stratification,fraction_edges_correctly_aligned,fraction_aligned_diffexp_edges,num_downstream_regulons,num_diffexp_regulons from cm_flows cmf join regulons r on cmf.regulon_id=r.id join mutations m on cmf.mutation_id=m.id join cm_flow_types cmft on cmf.cmf_type_id=cmft.id left outer join cmf_pathways as pw on pw.id=cmf.cmf_pathway_id left outer join genes mg on cmf.mutation_gene_id=mg.id join genes tf on cmf.regulator_id=tf.id where r.name=%s""", [regulon])
        return _make_causalflow_results(cursor)
    except:
        traceback.print_exc()
    finally:
        cursor.close()
        conn.close()


@app.route('/causalflows_for_regulator/<regulator>')
def causalflows_for_regulator(regulator):
    """all causal flows for the specified regulator"""
    conn = dbconn()
    cursor = conn.cursor(buffered=True)
    try:
        cursor.execute("""select cmf_id,r.name,m.name as mutation,cmft.name as cmftype,pw.name as pathway,mg.ensembl_id as mut_ensembl_id,mg.preferred as mut_preferred,tf.ensembl_id as tf_ensembl_id,tf.preferred as tf_preferred,regulon_mutation_regulator_role_id,regulon_mutation_regulator_pvalue,regulon_regulator_spearman_r,regulon_regulator_spearman_pvalue,regulon_t_statistic,regulon_log10_p_stratification,fraction_edges_correctly_aligned,fraction_aligned_diffexp_edges,num_downstream_regulons,num_diffexp_regulons from cm_flows cmf join regulons r on cmf.regulon_id=r.id join mutations m on cmf.mutation_id=m.id join cm_flow_types cmft on cmf.cmf_type_id=cmft.id left outer join cmf_pathways as pw on pw.id=cmf.cmf_pathway_id left outer join genes mg on cmf.mutation_gene_id=mg.id join genes tf on cmf.regulator_id=tf.id where tf.ensembl_id=%s""", [regulator])
        return _make_causalflow_results(cursor)
    except:
        traceback.print_exc()
    finally:
        cursor.close()
        conn.close()


@app.route('/causalflows_for_program/<program>')
def causalflows_for_program(program):
    """all causal flows for all regulons with the specified program"""
    conn = dbconn()
    cursor = conn.cursor(buffered=True)
    try:
        cursor.execute("""select cmf_id,r.name,m.name as mutation,cmft.name as cmftype,pw.name as pathway,mg.ensembl_id as mut_ensembl_id,mg.preferred as mut_preferred,tf.ensembl_id as tf_ensembl_id,tf.preferred as tf_preferred,regulon_mutation_regulator_role_id,regulon_mutation_regulator_pvalue,regulon_regulator_spearman_r,regulon_regulator_spearman_pvalue,regulon_t_statistic,regulon_log10_p_stratification,fraction_edges_correctly_aligned,fraction_aligned_diffexp_edges,num_downstream_regulons,num_diffexp_regulons from cm_flows cmf join regulons r on cmf.regulon_id=r.id join mutations m on cmf.mutation_id=m.id join cm_flow_types cmft on cmf.cmf_type_id=cmft.id left outer join cmf_pathways as pw on pw.id=cmf.cmf_pathway_id left outer join genes mg on cmf.mutation_gene_id=mg.id join genes tf on cmf.regulator_id=tf.id where r.id in (select regulon_id from regulon_programs rp join trans_programs p on rp.program_id=p.id where p.name=%s)""", [program])
        return _make_causalflow_results(cursor)
    except:
        traceback.print_exc()
    finally:
        cursor.close()
        conn.close()


@app.route('/causalflows_for_mutation/<mutation>')
def causalflows_for_mutation(mutation):
    """all causal flows for the specified mutation"""
    conn = dbconn()
    cursor = conn.cursor(buffered=True)
    try:
        cursor.execute("""select cmf_id,r.name,m.name as mutation,cmft.name as cmftype,pw.name as pathway,mg.ensembl_id as mut_ensembl_id,mg.preferred as mut_preferred,tf.ensembl_id as tf_ensembl_id,tf.preferred as tf_preferred,regulon_mutation_regulator_role_id,regulon_mutation_regulator_pvalue,regulon_regulator_spearman_r,regulon_regulator_spearman_pvalue,regulon_t_statistic,regulon_log10_p_stratification,fraction_edges_correctly_aligned,fraction_aligned_diffexp_edges,num_downstream_regulons,num_diffexp_regulons from cm_flows cmf join regulons r on cmf.regulon_id=r.id join mutations m on cmf.mutation_id=m.id join cm_flow_types cmft on cmf.cmf_type_id=cmft.id left outer join cmf_pathways as pw on pw.id=cmf.cmf_pathway_id left outer join genes mg on cmf.mutation_gene_id=mg.id join genes tf on cmf.regulator_id=tf.id where m.name=%s""", [mutation])
        return _make_causalflow_results(cursor)
    except:
        traceback.print_exc()
    finally:
        cursor.close()
        conn.close()


@app.route('/cfsearch/<term>')
def causal_flow_search(term):
    return jsonify(found="yes")


def _program_completions(cursor, prefix):
    cursor.execute('select distinct name from trans_programs where name like \"' + prefix + '%\" order by name')
    return [{'id': row[0], 'label': row[0], 'value': row[0]} for row in cursor.fetchall()]

@app.route('/completions/<term>')
def completions(term):
    """this is a function to serve the jquery autocomplete box"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        if term.startswith('P-'):
            completions = _program_completions(cursor, term)
        else:
            cursor.execute('select preferred from genes where preferred like %s', ["%s%%" % term])
            preferred = [{"id": row[0], "label": row[0], "value": row[0]} for row in cursor.fetchall()]
            cursor.execute('select ensembl_id from genes where ensembl_id like %s', ["%s%%" % term])
            ensembl = [{"id": row[0], "label": row[0], "value": row[0]} for row in cursor.fetchall()]

            completions = preferred + ensembl
        return jsonify(completions=completions)
    finally:
        cursor.close()
        conn.close()


@app.route('/mutation/<mutation_name>')
def mutation(mutation_name):
    """information for the specified mutation"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute("""select tfs.ensembl_id,tfs.preferred,bc.name,bmt.role_id,bc.cox_hazard_ratio
from regulon_mutation_regulator bmt join regulons bc on bmt.regulon_id=bc.id
join mutations m on m.id=bmt.mutation_id
join genes as tfs on tfs.id=bmt.regulator_id
where m.name=%s""",
                       [mutation_name])
        result = [{"regulator": tf, "regulator_preferred": tf_preferred if tf_preferred is not None else tf, "bicluster": bc,
                   "role": MUTATION_REGULATOR_ROLES[role],
                   "bc_cox_hazard_ratio": bc_cox_hazard_ratio,
                   "trans_program": ""}  # trans_program is N:M TODO
                  for tf, tf_preferred, bc, role,bc_cox_hazard_ratio in cursor.fetchall()]
        return jsonify(mutation=mutation_name, entries=result)
    finally:
        cursor.close()
        conn.close()


@app.route('/regulator/<tf_name>')
def regulator(tf_name):
    """information for the specified mutation"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute("""select bc.name,bt.role_id,bc.cox_hazard_ratio,mut.name
from regulon_regulator bt join regulons bc on bt.regulon_id=bc.id join genes as tfs on bt.regulator_id=tfs.id
join regulon_mutation_regulator bmt on bmt.regulon_id=bt.regulon_id and bmt.regulator_id=bt.regulator_id
join mutations mut on mut.id=bmt.mutation_id where tfs.ensembl_id=%s""",
                       [tf_name])
        result = [{
            "regulon": bc, "role": REGULATOR_REGULON_ROLES[role],
            "hazard_ratio": bc_hazard_ratio,
            "mutation": mut,
        } for bc, role, bc_hazard_ratio, mut in cursor.fetchall()]
        cursor.execute('select preferred from genes where ensembl_id=%s', [tf_name])
        row = cursor.fetchone()
        if row is not None and row[0] is not None:
            reg_preferred = row[0]
        else:
            reg_preferred = tf_name
        return jsonify(regulator=tf_name, regulator_preferred=reg_preferred,
                       entries=result)
    finally:
        cursor.close()
        conn.close()


@app.route('/program/<progname>')
def program(progname):
    """information for the specified program"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute("""
        select bc.id,name,bc.cox_hazard_ratio,num_genes,num_causal_flows
          from regulons bc
          join (select regulon_id, count(gene_id) as num_genes from regulon_genes
          group by regulon_id) as bcg on bc.id=bcg.regulon_id
          join (select regulon_id, count(*) as num_causal_flows
          from regulon_mutation_regulator
          group by regulon_id) as bccf on bc.id=bccf.regulon_id
          where bc.id in (select regulon_id from regulon_programs bp join trans_programs p on bp.program_id=p.id where p.name=%s)""",
                       [progname])

        regulons = []
        genes = []
        for regulon_id, name, cox_hazard_ratio, num_genes, num_causal_flows in cursor.fetchall():
            regulons.append({
                'regulon': regulon_id,
                'name': name,
                'cox_hazard_ratio': cox_hazard_ratio,
                'num_genes': num_genes,
                'num_causal_flows': num_causal_flows})
        cursor.execute("""
        select distinct ensembl_id,entrez_id,preferred from genes g join
        regulon_genes bg on g.id=bg.gene_id
        where bg.regulon_id in (select regulon_id from regulon_programs bp join trans_programs p on bp.program_id=p.id where p.name=%s)""", [progname])
        for ensembl_id, entrez_id, preferred in cursor.fetchall():
            genes.append({'ensembl_id': ensembl_id, 'entrez_id': entrez_id, 'preferred': preferred})
        return jsonify(regulons=regulons, genes=genes, num_genes=len(genes), num_regulons=len(regulons))
    finally:
        cursor.close()
        conn.close()

@app.route('/summary')
def summary():
    """model summary"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute('select count(*) from regulons')
        num_biclusters = cursor.fetchone()[0]
        cursor.execute('select count(*) from genes where is_regulator=1')
        num_regulators = cursor.fetchone()[0]
        cursor.execute('select count(*) from mutations')
        num_mutations = cursor.fetchone()[0]
        cursor.execute('select count(*) from cm_flows')
        num_causal_flows = cursor.fetchone()[0]
        cursor.execute('select count(*) from trans_programs')
        num_trans_programs = cursor.fetchone()[0]
        return jsonify(num_biclusters=num_biclusters,
                       num_regulators=num_regulators,
                       num_mutations=num_mutations,
                       num_causal_flows=num_causal_flows,
                       num_trans_programs=num_trans_programs)
    finally:
        cursor.close()
        conn.close()


@app.route('/gene_info/<gene>')
def gene_info(gene):
    """return all biclusters that contain this gene"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute('select entrez_id,ensembl_id,preferred from genes where entrez_id=%s or ensembl_id=%s or preferred=%s',
                       [gene, gene, gene])
        results = [(entrez_id, ensembl_id, preferred)
                   for entrez_id, ensembl_id, preferred in cursor.fetchall()]
        if len(results) > 0:
            ensembl_id = results[0][1]
            r = requests.get('https://rest.ensembl.org/lookup/id/%s?content-type=application/json;expand=1' % ensembl_id)
            if r.status_code == 200:
                ensdata = r.json()
                print(ensdata['description'])
            r = requests.get('https://rest.ensembl.org/xrefs/id/' + ensembl_id +
                             '?content-type=application/json;external_db=uniprot%;all_levels=1')
            if r.status_code == 200:
                xrefdata = r.json()
                uniprot_ids = [entry['primary_id'] for entry in xrefdata]
                if len(uniprot_ids) > 0:
                    uniprot_id = uniprot_ids[0]
                else:
                    uniprot_id = '-'
            if uniprot_id != '-':
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
            return jsonify(entrez_id=results[0][0], ensembl_id=ensembl_id,
                           preferred=results[0][2], description=ensdata['description'],
                           uniprot_id=uniprot_id, function=function)
        else:
            return jsonify(entrez_id="NA", ensembl_id="NA", preferred="NA", description='NA',
                           uniprot_id='NA', function='NA')
    finally:
        cursor.close()
        conn.close()


@app.route('/biclusters_for_gene/<gene_name>')
def biclusters_for_gene(gene_name):
    """return all biclusters that contain this gene"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute('select bc.name,cox_hazard_ratio from bicluster_genes bg join genes g on bg.gene_id=g.id join biclusters bc on bc.id=bg.bicluster_id where g.entrez_id=%s or g.ensembl_id=%s or g.preferred=%s',
                       [gene_name, gene_name, gene_name])
        return jsonify(biclusters=[{'cluster_id': row[0],
                                    'hazard_ratio': row[1]} for row in cursor.fetchall()])
    finally:
        cursor.close()
        conn.close()


@app.route('/bicluster_network/<cluster_id>')
def bicluster_network(cluster_id):
    """all genes in the specified bicluster"""
    conn = dbconn()
    cursor = conn.cursor()
    elements = [{"data": {"id": cluster_id}, "classes": "bicluster"}]
    try:
        edge_count = 0
        # retrieve the mutation and transcription factor nodes/edges
        # transcription factor -> bicluster
        cursor.execute('select tfs.ensembl_id,tfs.preferred,role from biclusters bc join regulon_regulator bt on bc.id=bt.regulon_id join genes tfs on tfs.id=bt.regulator_id where bc.name=%s', [cluster_id])
        for tf, tf_preferred, role in cursor.fetchall():
            tf = tf_preferred if tf_preferred is not None else tf
            elements.append({"data": {"id": tf}, "classes": "tf"})
            elements.append({"data": {"id": str(edge_count), "source": tf, "target": cluster_id},
                                      "classes": REGULATOR_REGULON_ROLES[role]})
            edge_count += 1

        # mutation role -> transcription factors
        cursor.execute('select m.name,tfs.ensembl_id,tfs.preferred,role from biclusters bc join bc_mutation_tf bmt on bc.id=bmt.bicluster_id join mutations m on m.id=bmt.mutation_id join genes tfs on tfs.id=bmt.tf_id where bc.name=%s', [cluster_id])
        for mutation, tf, tf_preferred, role in cursor.fetchall():
            tf = tf_preferred if tf_preferred is not None else tf
            elements.append({"data": {"id": mutation}, "classes": "mutation"})
            elements.append({"data": {"id": str(edge_count), "source": mutation, "target": tf},
                                      "classes": MUTATION_REGULATOR_ROLES[role].replace('-', '_') })
            edge_count += 1


        return jsonify(elements=elements)
    finally:
        cursor.close()
        conn.close()


@app.route('/bicluster_expressions/<cluster_id>')
def bicluster_expression_data(cluster_id):
    """returns data plot data in Highcharts format for bicluster expressions"""
    """this is actually box plot data, so the series needs to be a list of
    six tuples [condition_id, min, lower quartile, mean, upper quartile, max]]
    """
    conn = dbconn()
    cursor = conn.cursor()
    data = []
    try:
        cursor.execute('select p.name,median,min_value,max_value,lower_quartile,upper_quartile from bicluster_boxplot_data bbd join patients p on bbd.patient_id=p.id join biclusters bc on bbd.bicluster_id=bc.id where bc.name=%s order by median',
                       [cluster_id])
        for patient, median, minval, maxval, lower_quart, upper_quart in cursor.fetchall():
            data.append([patient, minval, lower_quart, median, upper_quart, maxval])
        return jsonify(data=data)
    finally:
        cursor.close()
        conn.close()


@app.route('/bicluster_enrichment/<cluster_id>')
def bicluster_enrichment(cluster_id):
    """returns barplot enrichment data for tumor subtypes in quitiles
    for the given bicluster"""
    conn = dbconn()
    cursor = conn.cursor()
    subtypes = ['g_cimp', 'proneural', 'neural', 'classical', 'mesenchymal', 'control']
    # series is gene -> list of values
    series  = defaultdict(list)
    # mockup some data for now (3 conditions)
    for s in subtypes:
        series[s] = [random.uniform(-10.0, 10.0) for i in range(5)]
    conds = ['All', 'All', 'All', 'All', 'All']
    return jsonify(expressions=series, conditions=conds)


MUTATION_ROLES = {1: 'down-regulates', 2: 'up-regulates'}
REGULATOR_ROLES = {1: 'activates', 2: 'represses'}

# TODO: This should actually be simply a pregenerated JSON file
@app.route('/causal_flow')
def causal_flow():
    """causal flow"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute("""select bc.name,mut.name,tfs.name,g.preferred,bmt.role,bc_tf.role_id,bc.cox_hazard_ratio,bgg.num_genes,bc.trans_program from bc_mutation_tf bmt join biclusters bc on bmt.bicluster_id=bc.id join mutations mut on bmt.mutation_id=mut.id join tfs on bmt.tf_id=tfs.id join bc_tf on bc.id=bc_tf.bicluster_id and tfs.id=bc_tf.tf_id join (select bc.id,count(bg.gene_id) as num_genes from biclusters bc join bicluster_genes bg on bc.id=bg.bicluster_id group by bc.id) as bgg on bc.id=bgg.id left join genes g on g.ensembl_id=tfs.name""")
        return jsonify(entries=[{
            'regulon': bc,
            'mutation': mut,
            'regulator': regulator,
            'regulator_preferred': tf_preferred if tf_preferred is not None else tf,
            'mutation_role': MUTATION_ROLES[mut_role],
            'regulator_role': REGULATOR_ROLES[tf_role],
            'hazard_ratio': hratio,
            'num_genes': ngenes,
            'trans_program': trans_program
        } for regulon,mut,regulator,tf_preferred,mut_role,tf_role,hratio,ngenes,trans_program in cursor.fetchall()])
    finally:
        cursor.close()
        conn.close()


"""
EXCELRA ACCESS FUNCTIONS
"""

@app.route('/diseases', methods=['POST'])
def diseases():
    reqdata = request.get_json()
    hr = reqdata['hr']
    regulon = reqdata['regulon'] if 'regulon' in reqdata else None
    mutation = reqdata['mutation'] if 'mutation' in reqdata else None
    regulator = reqdata['regulator'] if 'regulator' in reqdata else None

    conn = dbconn()
    cursor = conn.cursor()
    if mutation is None or mutation == 'All':
        canc_query = 'select distinct(c.name) from exc_cancer_mutation cm join exc_cancers c on cm.cancer_id=c.id where cm.hr=%s'
        cursor.execute(canc_query, [hr])
    else:
        canc_query = 'select distinct(c.name) from exc_cancer_mutation cm join exc_cancers c on cm.cancer_id=c.id join exc_mutations m on cm.mutation_id=m.id where cm.hr=%s and m.name=%s'
        cursor.execute(canc_query, [hr, mutation])
    cancers = [row[0] for row in cursor.fetchall()]

    # diseases are used in exc_disease_mutation, exc_disease_regulator, exc_disease_regulon
    if mutation is None or mutation == 'All':
        dis_query1 = 'select distinct(d.name) from exc_disease_mutation dm join exc_diseases d on dm.disease_id=d.id where dm.hr=%s'
        cursor.execute(dis_query1, [hr])
    else:
        dis_query1 = 'select distinct(d.name) from exc_disease_mutation dm join exc_diseases d on dm.disease_id=d.id join exc_mutations m on dm.mutation_id=m.id where dm.hr=%s and m.name=%s'
        cursor.execute(dis_query1, [hr, mutation])
    diseases1 = [row[0] for row in cursor.fetchall()]

    if regulator is None or regulator == 'All':
        dis_query2 = 'select distinct(d.name) from exc_disease_regulator dr join exc_diseases d on dr.disease_id=d.id where dr.hr=%s'
        cursor.execute(dis_query2, [hr])
    else:
        dis_query2 = 'select distinct(d.name) from exc_disease_regulator dr join exc_diseases d on dr.disease_id=d.id join exc_regulators r on dr.regulator_id=r.id where dr.hr=%s and r.name=%s'
        cursor.execute(dis_query2, [hr, regulator])
    diseases2 = [row[0] for row in cursor.fetchall()]

    if regulon is None:
        dis_query3 = 'select distinct(d.name) from exc_disease_regulon dr join exc_diseases d on dr.disease_id=d.id where dr.hr=%s'
        cursor.execute(dis_query3, [hr])
    else:
        dis_query3 = 'select distinct(d.name) from exc_disease_regulon dr join exc_diseases d on dr.disease_id=d.id join exc_regulons r on dr.regulon_id=r.id where dr.hr=%s and r.name=%s'
        cursor.execute(dis_query3, [hr, regulon])
    diseases3 = [row[0] for row in cursor.fetchall()]

    disease_cancers = set(cancers + diseases1 + diseases2 + diseases3)
    return jsonify(status='ok', diseases=sorted(disease_cancers))


@app.route('/mutations', methods=['POST'])
def mutations():
    reqdata = request.get_json()
    hr = reqdata['hr']
    # mutations are dependent on drugs, diseases and regulators
    drug = reqdata['drug'] if 'drug' in reqdata else None
    disease = reqdata['disease'] if 'disease' in reqdata else None
    regulator = reqdata['regulator'] if 'regulator' in reqdata else None
    conn = dbconn()
    cursor = conn.cursor()

    if disease is None or disease == 'All myelomas' or disease == 'All Cancers':
        cursor.execute('select m.name from exc_disease_mutation dm join exc_mutations m on dm.mutation_id=m.id where dm.hr=%s', [hr])
    else:
        cursor.execute('select m.name from exc_disease_mutation dm join exc_mutations m on dm.mutation_id=m.id join exc_diseases d on dm.disease_id=d.id where dm.hr=%s and d.name=%s',
                       [hr, disease])

    mutations1 = [row[0] for row in cursor.fetchall()]

    if drug is not None:
        cursor.execute('select m.name from exc_mutation_drug md join exc_mutations m on md.mutation_id=m.id join exc_drugs d on md.drug_id=d.id where md.hr=%s and d.name=%s', [hr, drug])
    else:
        cursor.execute('select m.name from exc_mutation_drug md join exc_mutations m on md.mutation_id=m.id where md.hr=%s', [hr])
    mutations2 = [row[0] for row in cursor.fetchall()]

    if regulator is None or regulator == 'All':
        cursor.execute('select m.name from exc_mutation_regulator mr join exc_mutations m on mr.mutation_id=m.id where mr.hr=%s', [hr])
    else:
        cursor.execute('select m.name from exc_mutation_regulator mr join exc_mutations m on mr.mutation_id=m.id join exc_regulators r on mr.regulator_id=r.id where mr.hr=%s and r.name=%s',
                       [hr, regulator])
    mutations3 = [row[0] for row in cursor.fetchall()]
    return jsonify(status='ok', mutations=sorted(set(mutations1 + mutations2 + mutations3)))


@app.route('/regulators', methods=['POST'])
def regulators():
    reqdata = request.get_json()
    hr = reqdata['hr']
    disease = reqdata['disease'] if 'disease' in reqdata else None
    mutation = reqdata['mutation'] if 'mutation' in reqdata else None

    # filter by the other attributes, e.g. disease,
    # in tables exc_disease_regulator, exc_mutation_regulator
    conn = dbconn()
    cursor = conn.cursor()
    if disease is None or disease == 'All myelomas' or disease == 'All Cancers':
        cursor.execute('select r.name from exc_disease_regulator dr join exc_regulators r on dr.regulator_id=r.id where dr.hr=%s', [hr])
    else:
        cursor.execute('select r.name from exc_disease_regulator dr join exc_regulators r on dr.regulator_id=r.id join exc_diseases d on dr.disease_id=d.id where dr.hr=%s and d.name=%s', [hr, disease])

    regulators1 = [row[0] for row in cursor.fetchall()]

    if mutation is None or mutation == 'All':
        cursor.execute('select r.name from exc_mutation_regulator mr join exc_regulators r on mr.regulator_id=r.id where mr.hr=%s', [hr])
    else:
        cursor.execute('select r.name from exc_mutation_regulator mr join exc_regulators r on mr.regulator_id=r.id join exc_mutations m on mr.mutation_id=m.id where mr.hr=%s and m.name=%s',
                       [hr, mutation])
    regulators2 = [row[0] for row in cursor.fetchall()]
    return jsonify(status='ok', regulators=sorted(set(regulators1 + regulators2)))


@app.route('/regulons', methods=['POST'])
def regulons():
    reqdata = request.get_json()
    hr = reqdata['hr']
    disease = reqdata['disease'] if 'disease' in reqdata else None

    # in tables exc_disease_regulon
    conn = dbconn()
    cursor = conn.cursor()
    if disease is None or disease == 'All myelomas' or disease == 'All Cancers':
        cursor.execute('select r.name from exc_disease_regulon dr join exc_regulons r on dr.regulon_id=r.id where dr.hr=%s order by name',
                   [hr])
    else:
        cursor.execute('select r.name from exc_disease_regulon dr join exc_regulons r on dr.regulon_id=r.id join exc_diseases d on dr.disease_id=d.id where dr.hr=%s and d.name=%s order by name',
                   [hr, disease])

    regulons = [row[0] for row in cursor.fetchall()]
    return jsonify(status='ok', regulons=regulons)


@app.route('/drugs', methods=['POST'])
def drugs():
    reqdata = request.get_json()
    hr = reqdata['hr']
    mutation = reqdata['mutation'] if 'mutation' in reqdata else None

    # in tables exc_mutation_drug
    conn = dbconn()
    cursor = conn.cursor()
    if mutation is None or mutation == 'All':
        cursor.execute('select d.name from exc_mutation_drug md join exc_drugs d on md.drug_id=d.id where md.hr=%s order by name',
                       [hr])
    else:
        cursor.execute('select d.name from exc_mutation_drug md join exc_drugs d on md.drug_id=d.id join exc_mutations m on md.mutation_id=m.id where md.hr=%s and m.name=%s order by name',
                       [hr, mutation])

    drugs = [row[0] for row in cursor.fetchall()]
    return jsonify(status='ok', drugs=drugs)


@app.route('/pmid_counts/<hr>/<disease>/<mutation>/<regulator>/<regulon>/<drug>')
def pmid_counts(hr, disease, mutation, regulator, regulon, drug):
    conn = dbconn()
    cursor = conn.cursor()
    return jsonify(status='ok',
                   num_disease_mutation_pmids=len(datasets.disease_mutation(cursor, hr, disease, mutation)),
                   num_disease_regulator_pmids=len(datasets.disease_regulator(cursor, hr, disease, regulator)),
                   num_disease_regulon_pmids=len(datasets.disease_regulon(cursor, hr, disease, regulon)),
                   num_mutation_regulator_pmids=len(datasets.mutation_regulator(cursor, hr, mutation, regulator)),
                   num_mutation_drug_pmids=len(datasets.mutation_drug(cursor, hr, mutation, drug)))


@app.route('/search_pmid_counts/<hr>', methods=['POST'])
def search_pmid_counts(hr):
    conn = dbconn()
    cursor = conn.cursor()
    reqdata = request.get_json()
    search_term = reqdata['search']
    print(search_term)
    return jsonify(status="ok",
                   num_disease_mutation_pmids=len(datasets.search_disease_mutation(cursor, hr, search_term)),
                   num_disease_regulator_pmids=len(datasets.search_disease_regulator(cursor, hr, search_term)),
                   num_disease_regulon_pmids=len(datasets.search_disease_regulon(cursor, hr, search_term)),
                   num_mutation_regulator_pmids=len(datasets.search_mutation_regulator(cursor, hr, search_term)),
                   num_mutation_drug_pmids=len(datasets.search_mutation_drug(cursor, hr, search_term)))


def fetch_articles(pmids):
    pmid_str = ','.join(map(str, pmids))
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=%s&retmode=xml" % pmid_str

    r = requests.get(url)
    try:
        root = ET.fromstring(r.content)
        result = []
        for item in root.findall('./PubmedArticle/MedlineCitation'):
            pmid = item.find('PMID').text
            article = item.find('Article')
            title = article.find('ArticleTitle').text
            try:
                abstract = article.find('Abstract/AbstractText').text
            except:
                abstract = ''
            pubdate_str = ''
            try:
                for pubdate in article.findall('.//PubDate'):
                    pubyear = pubdate.find('Year').text
                    pubmonth = pubdate.find('Month').text
                    pubdate_str = '%s/%s' % (pubyear, pubmonth)
            except:
                pass

            result.append({'pmid': pmid, 'title': title, 'abstract': abstract,
                           'pubdate': pubdate_str})

        return result
    except:
        raise


def batch_results(payload, pmids):
    items_per_page = payload['itemsPerPage']
    page = payload['page']

    result = []
    if len(pmids) > 0:
        total = len(pmids)
        if len(pmids) > items_per_page:
            offset = (page - 1) * items_per_page
            pmids = pmids[offset:offset + items_per_page]
        return fetch_articles(pmids)
    return result


def markup_assocs(content, assocs):
    for assoc in assocs:
        regexp = re.compile(re.escape(assoc), re.IGNORECASE)
        content = regexp.sub('<span class="marked">' + assoc + '</span>', content)
    return content

@app.route('/cancer_mutation_docs/<hr>/<cancer>/<mutation>', methods=['POST'])
def cancer_mutation_docs(hr, cancer, mutation):
    conn = dbconn()
    cursor = conn.cursor()
    pmids = datasets.cancer_mutation(cursor, hr, cancer, mutation)
    total = len(pmids)
    result = batch_results(request.get_json(), pmids)
    for entry in result:
        assocs = []
        if cancer != 'All':
            assocs.append(cancer)
        if mutation != 'All':
            assocs.append(mutation)

        entry['abstract'] = markup_assocs(entry['abstract'], assocs)
        entry['assocs'] = '->'.join(list(assocs))
    return jsonify(status='ok', total=total, data=result)


def __make_search_results(pmids, search_term):
    total = len(pmids)
    result = batch_results(request.get_json(), pmids)
    assocs = [search_term]
    for entry in result:
        entry['assocs'] = search_term
        entry['abstract'] = markup_assocs(entry['abstract'], assocs)
    return result, total


@app.route('/disease_mutation_docs/<hr>/<disease>/<mutation>', methods=['POST'])
def disease_mutation_docs(hr, disease, mutation):
    conn = dbconn()
    cursor = conn.cursor()
    pmids = datasets.disease_mutation(cursor, hr, disease, mutation)
    total = len(pmids)
    result = batch_results(request.get_json(), pmids)
    for entry in result:
        assocs = []
        if disease != 'All':
            assocs.append(disease)
        if mutation != 'All':
            assocs.append(mutation)

        entry['assocs'] = '->'.join(list(assocs))
        entry['abstract'] = markup_assocs(entry['abstract'], assocs)
    return jsonify(status='ok', total=total, data=result)


@app.route('/disease_mutation_search/<hr>/<search_term>', methods=['POST'])
def disease_mutation_search(hr, search_term):
    conn = dbconn()
    cursor = conn.cursor()
    pmids = datasets.search_disease_mutation(cursor, hr, search_term)
    result, total = __make_search_results(pmids, search_term)
    return jsonify(status='ok', total=total, data=result)


@app.route('/disease_regulator_docs/<hr>/<disease>/<regulator>', methods=['POST'])
def disease_regulator_docs(hr, disease, regulator):
    conn = dbconn()
    cursor = conn.cursor()
    pmids = datasets.disease_regulator(cursor, hr, disease, regulator)
    total = len(pmids)
    result = batch_results(request.get_json(), pmids)
    for entry in result:
        assocs = []
        if disease != 'All':
            assocs.append(disease)
        if regulator != 'All':
            assocs.append(regulator)

        entry['assocs'] = '->'.join(list(assocs))
        try:
            entry['abstract'] = markup_assocs(entry['abstract'], assocs)
        except:
            traceback.print_exc()
    return jsonify(status='ok', total=total, data=result)


@app.route('/disease_regulator_search/<hr>/<search_term>', methods=['POST'])
def disease_regulator_search(hr, search_term):
    conn = dbconn()
    cursor = conn.cursor()
    pmids = datasets.search_disease_regulator(cursor, hr, search_term)
    result, total = __make_search_results(pmids, search_term)
    return jsonify(status='ok', total=total, data=result)


@app.route('/disease_regulon_docs/<hr>/<disease>/<regulon>', methods=['POST'])
def disease_regulon_docs(hr, disease, regulon):
    conn = dbconn()
    cursor = conn.cursor()
    pmids = datasets.disease_regulon(cursor, hr, disease, regulon)
    total = len(pmids)
    result = batch_results(request.get_json(), pmids)
    for entry in result:
        assocs = []
        if disease != 'All':
            assocs.append(disease)
        if regulon != 'All':
            assocs.append(regulon)

        entry['assocs'] = '->'.join(list(assocs))
        entry['abstract'] = markup_assocs(entry['abstract'], assocs)
    return jsonify(status='ok', total=total, data=result)


@app.route('/disease_regulon_search/<hr>/<search_term>', methods=['POST'])
def disease_regulon_search(hr, search_term):
    conn = dbconn()
    cursor = conn.cursor()
    pmids = datasets.search_disease_regulon(cursor, hr, search_term)
    result, total = __make_search_results(pmids, search_term)
    return jsonify(status='ok', total=total, data=result)


@app.route('/mutation_regulator_docs/<hr>/<mutation>/<regulator>', methods=['POST'])
def mutation_regulator_docs(hr, mutation, regulator):
    conn = dbconn()
    cursor = conn.cursor()
    pmids = datasets.mutation_regulator(cursor, hr, mutation, regulator)
    total = len(pmids)
    result = batch_results(request.get_json(), pmids)
    for entry in result:
        assocs = []
        if mutation != 'All':
            assocs.append(mutation)
        if regulator != 'All':
            assocs.append(regulator)

        entry['assocs'] = '->'.join(list(assocs))
        entry['abstract'] = markup_assocs(entry['abstract'], assocs)
    return jsonify(status='ok', total=total, data=result)


@app.route('/mutation_regulator_search/<hr>/<search_term>', methods=['POST'])
def mutation_regulator_search(hr, search_term):
    conn = dbconn()
    cursor = conn.cursor()
    pmids = datasets.search_mutation_regulator(cursor, hr, search_term)
    result, total = __make_search_results(pmids, search_term)
    return jsonify(status='ok', total=total, data=result)


@app.route('/mutation_drug_docs/<hr>/<mutation>/<drug>', methods=['POST'])
def mutation_drug_docs(hr, mutation, drug):
    conn = dbconn()
    cursor = conn.cursor()
    pmids = datasets.mutation_drug(cursor, hr, mutation, drug)
    total = len(pmids)
    result = batch_results(request.get_json(), pmids)
    for entry in result:
        assocs = []
        if mutation != 'All':
            assocs.append(mutation)
        if drug != 'All':
            assocs.append(drug)

        entry['assocs'] = '->'.join(list(assocs))
        entry['abstract'] = markup_assocs(entry['abstract'], assocs)
    return jsonify(status='ok', total=total, data=result)



@app.route('/mutation_drug_search/<hr>/<search_term>', methods=['POST'])
def mutation_drug_search(hr, search_term):
    conn = dbconn()
    cursor = conn.cursor()
    pmids = datasets.search_mutation_drug(cursor, hr, search_term)
    result, total = __make_search_results(pmids, search_term)
    return jsonify(status='ok', total=total, data=result)


if __name__ == '__main__':
    app.debug = True
    app.secret_key = 'trstrestnorgp654g'
    app.run(host='0.0.0.0', debug=True)
