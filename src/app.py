#!/usr/bin/env python

import logging
import json
import os
import traceback
from collections import defaultdict
import mysql.connector

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

HALLMARK_IMGS = {
    "Self_Sufficiency_In_Growth_Signals": "self_sufficiency_in_growth_signals.gif",
    "Insensitivity_To_Antigrowth_Signals": "insensitivity_to_antigrowth_signals.gif",
    "Evading_Apoptosis": "evading_apoptosis.gif",
    "Limitless_Replicative_Potential": "limitless_replicative_potential.gif",
    "Sustained_Angiogenesis": "sustained_angiogenesis.gif",
    "Tissue_Invasion_And_Metastasis": "tissue_invasion_and_metastasis.gif",
    "Genome_Instability_And_Mutation": "genome_instability_and_mutation.gif",
    "Tumor_Promoting_Inflammation": "tumor_promoting_inflammation.gif",
    "Reprogramming_Energy_Metabolism": "reprogramming_energy_metabolism.gif",
    "Evading_ImmuneDetection": "evading_immune_detection.gif"
}

HALLMARK_NAMES = {
    "Self_Sufficiency_In_Growth_Signals": "Self sufficiency in growth signal",
    "Insensitivity_To_Antigrowth_Signals": "Insensitivity to antigrowth signal",
    "Evading_Apoptosis": "Evading apoptosis",
    "Limitless_Replicative_Potential": "Limitless replicative potential",
    "Sustained_Angiogenesis": "Sustained angiogenesis",
    "Tissue_Invasion_And_Metastasis": "Tissue invasion and metastasis",
    "Genome_Instability_And_Mutation": "Genome instability and mutation",
    "Tumor_Promoting_Inflammation": "Tumor-promoting inflammation",
    "Reprogramming_Energy_Metabolism": "Reprogramming energy metabolism",
    "Evading_ImmuneDetection": "Evading immune detection"
}

@app.route('/regulon/<regulon>')
def regulon_info(regulon):
    conn = dbconn()
    cursor = conn.cursor(buffered=True)
    try:
        cursor.execute("select id,cox_hazard_ratio from regulons where name=%s", [regulon])
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
        cursor.execute("""select g.entrez_id,g.ensembl_id,g.preferred,g.uniprot_id,g.uniprot_function,g.ens_description
from regulons bc join regulon_genes bcg on bc.id=bcg.regulon_id
join genes g on g.id=bcg.gene_id where bc.name=%s""",
                       [regulon])

        genes = []
        for entrez_id, ensembl_id, preferred, uniprot_id, function, description in cursor.fetchall():
            genes.append({
                "entrez_id": entrez_id,
                "ensembl_id": ensembl_id,
                "preferred": preferred,
                "uniprot_id": uniprot_id,
                "function": function,
                "description": description
            })

        # Regulon hallmarks
        cursor.execute("""select distinct h.name from regulon_programs rp join program_hallmarks ph on rp.program_id=ph.program_id join hallmarks h on h.id=ph.hallmark_id join regulons r on rp.regulon_id=rp.regulon_id where r.name=%s and ph.wang_score >= 0.8""", [regulon])
        hallmarks = [{
            "name": HALLMARK_NAMES[row[0]],
            "image": HALLMARK_IMGS[row[0]]
        } for row in cursor.fetchall()]

        # transcriptional programs
        cursor.execute("""select tp.name
from regulons reg join regulon_programs rp on reg.id=rp.regulon_id
join trans_programs tp on tp.id=rp.program_id where reg.name=%s""", [regulon])
        programs = [row[0] for row in cursor.fetchall()]

        return jsonify(regulon=regulon,
                       hazard_ratio=hazard_ratio,
                       regulon_regulators=regulon_regulators,
                       genes=genes,
                       hallmarks=hallmarks,
                       program=programs[0],  # for now, just return 1 element, it's actually 1:1
                       num_causal_flows=num_causalflows)
    except:
        traceback.print_exc()
    finally:
        cursor.close()
        conn.close()

def _make_causalflow_row(conn, row):
    (cmf_id, regulon, mutation, cmf_type, pathway,
     mutgene_ensembl, mutgene_symbol,
     regulator_ensembl, regulator_symbol, mutation_regulator_role,
     regulon_mutation_regulator_pvalue,
     regulon_regulator_spearman_r, regulon_regulator_spearman_pvalue,
     regulon_t_statistic, regulon_log10_p_stratification,
     fraction_edges_correctly_aligned, fraction_aligned_diffexp_edges,
     num_downstream_regulons, num_diffexp_regulons) = row

    # Addon: extract drugs associated with the row according to rules
    regulator_drugs = []
    regulon_drugs = []
    cur = conn.cursor()
    try:
        # program related drugs
        #cur.execute("select distinct d.name from drugs d join drug_programs dp on d.id=dp.drug_id join trans_programs p on p.id=dp.program_id join target_type_models ttm on d.target_type_model_id=ttm.id where ttm.name='Regulon_Regulator' and dp.program_id in (select rp.program_id from regulon_programs rp join regulons r on rp.regulon_id=r.id where r.name=%s)", [regulon])
        cur.execute("select distinct d.name from drugs d join drug_regulons dr on d.id=dr.drug_id join regulons r on r.id=dr.regulon_id join target_type_models ttm on d.target_type_model_id=ttm.id where ttm.name='Regulon_Regulator' and r.name=%s", [regulon])

        for row in cur.fetchall():
            regulator_drugs.append(row[0])

        regulon_drugs = _query_regulon_drugs(cur, regulon)

        #cur.execute("select distinct d.name from drugs d join drug_regulons dr on d.id=dr.drug_id join regulons r on r.id=dr.regulon_id join target_type_models ttm on d.target_type_model_id=ttm.id where ttm.name='Regulon_Gene'")
    finally:
        cur.close()

    return {
        "cmf_id": cmf_id,
        "regulon": regulon,
        "cmf_type": cmf_type,
        "pathway": pathway,
        "mutation_gene_ensembl": mutgene_ensembl,
        "mutation_gene_symbol": mutgene_symbol,
        "mutation": mutation,
        "regulator": regulator_ensembl,
        "regulator_preferred": regulator_symbol if regulator_symbol is not None else regulator_ensembl,
        "regulator_role": REGULATOR_REGULON_ROLES[1 if regulon_regulator_spearman_r < 0 else 2],
        "mutation_role": MUTATION_REGULATOR_ROLES[mutation_regulator_role],
        "regulator_pvalue": regulon_mutation_regulator_pvalue,
        "regulator_spearman_r": regulon_regulator_spearman_r,
        "regulator_spearman_pvalue": regulon_regulator_spearman_pvalue,
        "regulon_t_statistic": regulon_t_statistic,
        "regulon_log10_p_stratification": regulon_log10_p_stratification,
        "fraction_edges_correctly_aligned": fraction_edges_correctly_aligned,
        "fraction_aligned_diffexp_edges": fraction_aligned_diffexp_edges,
        "num_downstream_regulons": num_downstream_regulons,
        "num_diffexp_regulons": num_diffexp_regulons,
        "regulator_drugs": regulator_drugs,
        "regulon_drugs": regulon_drugs
    }


def _query_regulon_drugs(cur, regulon):
    result = []
    #cur.execute("select distinct d.name from drugs d join drug_targets dt on d.id=dt.drug_id join genes g on g.id=dt.gene_id join target_type_models ttm on d.target_type_model_id=ttm.id where ttm.name='Regulon_Gene' and g.id in (select g.id from genes g join regulon_genes rg on g.id=rg.gene_id join regulons r on r.id=rg.regulon_id where r.name=%s)", [regulon])
    cur.execute("select g.id from genes g join regulon_genes rg on g.id=rg.gene_id join regulons r on r.id=rg.regulon_id where r.name=%s", [regulon])
    gene_ids = [str(row[0]) for row in cur.fetchall()]
    if len(gene_ids) > 0:
        gids = "(%s)" % ','.join(gene_ids)
        query = """select distinct d.name,d.approved_symbol,d.approved,d.max_trial_phase,d.max_phase_gbm,dtp.name,
at.name,moa.name
from drugs d join drug_targets dt on d.id=dt.drug_id join genes g on g.id=dt.gene_id
join target_type_models ttm on d.target_type_model_id=ttm.id
join drug_types dtp on dtp.id=d.drug_type_id
join action_types at on at.id=d.action_type_id
join mechanisms_of_action moa on moa.id=d.mechanism_of_action_id
where ttm.name='Regulon_Gene' and g.id in """ + gids
        cur.execute(query)
        for (name, approved_symbol, approved, max_trial_phase,
             max_gbm_phase, drug_type, action_type, mechanism_of_action) in cur.fetchall():
            result.append({
                "name": name,
                "approved_symbol": approved_symbol,
                "max_trial_phase": max_trial_phase if max_trial_phase is not None else 0,
                "max_gbm_phase": max_gbm_phase if max_gbm_phase is not None else 0,
                "drug_type": drug_type,
                "action_type": action_type,
                "mechanism_of_action": mechanism_of_action
            })
    return result

@app.route('/regulon_drugs/<regulon>')
def regulon_drugs(regulon):
    """all drugs for the specified regulon"""
    conn = dbconn()
    cursor = conn.cursor(buffered=True)
    try:
        drugs = _query_regulon_drugs(cursor, regulon)
        return jsonify(drugs=drugs)
    finally:
        cursor.close()
        conn.close()


def _make_causalflow_results(conn, cursor):
    cm_flows = [_make_causalflow_row(conn, row) for row in cursor.fetchall()]
    return jsonify(cm_flows=cm_flows)


@app.route('/causalflows_for_regulon/<regulon>')
def causalflows_for_regulon(regulon):
    """all causal flows for the specified regulon"""
    conn = dbconn()
    cursor = conn.cursor(buffered=True)
    try:
        cursor.execute("""select cmf_id,r.name,m.name as mutation,cmft.name as cmftype,pw.name as pathway,mg.ensembl_id as mut_ensembl_id,mg.preferred as mut_preferred,tf.ensembl_id as tf_ensembl_id,tf.preferred as tf_preferred,regulon_mutation_regulator_role_id,regulon_mutation_regulator_pvalue,regulon_regulator_spearman_r,regulon_regulator_spearman_pvalue,regulon_t_statistic,regulon_log10_p_stratification,fraction_edges_correctly_aligned,fraction_aligned_diffexp_edges,num_downstream_regulons,num_diffexp_regulons from cm_flows cmf join regulons r on cmf.regulon_id=r.id join mutations m on cmf.mutation_id=m.id join cm_flow_types cmft on cmf.cmf_type_id=cmft.id left outer join cmf_pathways as pw on pw.id=cmf.cmf_pathway_id left outer join genes mg on cmf.mutation_gene_id=mg.id join genes tf on cmf.regulator_id=tf.id where r.name=%s""", [regulon])
        return _make_causalflow_results(conn, cursor)
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
        cursor.execute("""select cmf_id,r.name,m.name as mutation,cmft.name as cmftype,pw.name as pathway,mg.ensembl_id as mut_ensembl_id,mg.preferred as mut_preferred,tf.ensembl_id as tf_ensembl_id,tf.preferred as tf_preferred,regulon_mutation_regulator_role_id,regulon_mutation_regulator_pvalue,regulon_regulator_spearman_r,regulon_regulator_spearman_pvalue,regulon_t_statistic,regulon_log10_p_stratification,fraction_edges_correctly_aligned,fraction_aligned_diffexp_edges,num_downstream_regulons,num_diffexp_regulons from cm_flows cmf join regulons r on cmf.regulon_id=r.id join mutations m on cmf.mutation_id=m.id join cm_flow_types cmft on cmf.cmf_type_id=cmft.id left outer join cmf_pathways as pw on pw.id=cmf.cmf_pathway_id left outer join genes mg on cmf.mutation_gene_id=mg.id join genes tf on cmf.regulator_id=tf.id where tf.ensembl_id=%s or tf.preferred=%s""", [regulator, regulator])
        return _make_causalflow_results(conn, cursor)
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
        return _make_causalflow_results(conn, cursor)
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
        return _make_causalflow_results(conn, cursor)
    except:
        traceback.print_exc()
    finally:
        cursor.close()
        conn.close()


@app.route('/causalflows_with_mutation_in/<gene>')
def causalflows_with_mutation_in(gene):
    """all causal flows for the specified mutation"""
    conn = dbconn()
    cursor = conn.cursor(buffered=True)
    try:
        cursor.execute("""select cmf_id,r.name,m.name as mutation,cmft.name as cmftype,pw.name as pathway,mg.ensembl_id as mut_ensembl_id,mg.preferred as mut_preferred,tf.ensembl_id as tf_ensembl_id,tf.preferred as tf_preferred,regulon_mutation_regulator_role_id,regulon_mutation_regulator_pvalue,regulon_regulator_spearman_r,regulon_regulator_spearman_pvalue,regulon_t_statistic,regulon_log10_p_stratification,fraction_edges_correctly_aligned,fraction_aligned_diffexp_edges,num_downstream_regulons,num_diffexp_regulons from cm_flows cmf join regulons r on cmf.regulon_id=r.id join mutations m on cmf.mutation_id=m.id join cm_flow_types cmft on cmf.cmf_type_id=cmft.id left outer join cmf_pathways as pw on pw.id=cmf.cmf_pathway_id left outer join genes mg on cmf.mutation_gene_id=mg.id join genes tf on cmf.regulator_id=tf.id where mg.preferred=%s or mg.ensembl_id=%s""", [gene, gene])
        return _make_causalflow_results(conn, cursor)
    except:
        traceback.print_exc()
    finally:
        cursor.close()
        conn.close()


@app.route('/causalflows_for_regulons_containing/<gene>')
def causalflows_for_regulons_containing(gene):
    """all causal flows for regulons that contain the specified gene"""
    conn = dbconn()
    cursor = conn.cursor(buffered=True)
    try:
        # NOTE: This was originally a single query, and it worked on the development
        # environment, but crashed with an unknown server error and NO log message.
        # Splitting this up into several queries generally should NOT be
        # done, but it's the workaround for now
        cursor.execute("select distinct r.id from regulons r join regulon_genes rg on r.id=rg.regulon_id join genes g on rg.gene_id=g.id where g.preferred=%s or g.ensembl_id=%s", [gene, gene])
        regulon_ids = [row[0] for row in cursor.fetchall()]
        cm_flows = []
        for regulon_id in regulon_ids:
            cursor.execute("select cmf_id,r.name,m.name as mutation,cmft.name as cmftype,pw.name as pathway,mg.ensembl_id as mut_ensembl_id,mg.preferred as mut_preferred,tf.ensembl_id as tf_ensembl_id,tf.preferred as tf_preferred,regulon_mutation_regulator_role_id,regulon_mutation_regulator_pvalue,regulon_regulator_spearman_r,regulon_regulator_spearman_pvalue,regulon_t_statistic,regulon_log10_p_stratification,fraction_edges_correctly_aligned,fraction_aligned_diffexp_edges,num_downstream_regulons,num_diffexp_regulons from cm_flows cmf join regulons r on cmf.regulon_id=r.id join mutations m on cmf.mutation_id=m.id join cm_flow_types cmft on cmf.cmf_type_id=cmft.id left outer join cmf_pathways as pw on pw.id=cmf.cmf_pathway_id left outer join genes mg on cmf.mutation_gene_id=mg.id join genes tf on cmf.regulator_id=tf.id where r.id=%s", [regulon_id])
            for row in cursor.fetchall():
                cm_flows.append(_make_causalflow_row(conn, row))
        return jsonify(cm_flows=cm_flows)
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


@app.route('/mutations')
def mutations():
    """list of mutations"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute("""select distinct m.ensembl_id,m.preferred,m.entrez_id,mut.name from genes as m join cm_flows cmf on cmf.mutation_gene_id=m.id join mutations mut on cmf.mutation_id=mut.id where m.is_mutation=1""")
        result = [{"ensembl": mut_ensembl if mut_ensembl is not None else '-',
                   "preferred": mut_preferred if mut_preferred is not None else '-',
                   "entrez": mut_entrez if mut_entrez is not None else '-',
                   "mutation": mutation}
                  for mut_ensembl, mut_preferred, mut_entrez, mutation in cursor.fetchall()]
        return jsonify(entries=result)
    finally:
        cursor.close()
        conn.close()


@app.route('/regulators')
def regulators():
    """list of regulators"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute("""select distinct r.ensembl_id,r.preferred,r.entrez_id from genes as r where r.is_regulator=1 order by r.preferred""")
        result = [{"ensembl": ensembl if ensembl is not None else '-',
                   "preferred": preferred if preferred is not None else '-',
                   "entrez": entrez if entrez is not None else '-'}
                  for ensembl, preferred, entrez in cursor.fetchall()]
        return jsonify(entries=result)
    finally:
        cursor.close()
        conn.close()


@app.route('/regulons')
def regulons():
    """list of regulons"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute("""select distinct r.name,p.name,count(distinct rg.gene_id),count(distinct rr.regulator_id),count(distinct cmf.id) from regulons r join regulon_programs rp on r.id=rp.regulon_id join trans_programs p on p.id=rp.program_id join regulon_genes rg on r.id=rg.regulon_id join regulon_regulator rr on rr.regulon_id=r.id join cm_flows cmf on cmf.regulon_id=r.id group by p.id,rg.regulon_id, rr.regulon_id,cmf.regulon_id order by r.name""")
        result = [{"regulon": regulon,
                   "program": program,
                   "num_genes": num_genes,
                   "num_regulators": num_regulators,
                   "num_causalflows": num_causalflows}
                  for regulon, program, num_genes, num_regulators, num_causalflows in cursor.fetchall()]
        return jsonify(entries=result)
    finally:
        cursor.close()
        conn.close()


@app.route('/gene_regulons/<gene>')
def gene_regulons(gene):
    """list of regulons for a gene"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute("""select distinct r.name,p.name,count(distinct rg.gene_id),count(distinct rr.regulator_id),count(distinct cmf.id)
        from regulons r join regulon_programs rp on r.id=rp.regulon_id
        join trans_programs p on p.id=rp.program_id
        join regulon_genes rg on r.id=rg.regulon_id
        join regulon_regulator rr on rr.regulon_id=r.id
        join cm_flows cmf on cmf.regulon_id=r.id
        join genes g on g.id=rg.gene_id
        where g.preferred=%s
        group by p.id,rg.regulon_id, rr.regulon_id,cmf.regulon_id order by r.name""", [gene])

        result = [{"regulon": regulon,
                   "program": program,
                   "num_genes": num_genes,
                   "num_regulators": num_regulators,
                   "num_causalflows": num_causalflows}
                  for regulon, program, num_genes, num_regulators, num_causalflows in cursor.fetchall()]
        return jsonify(entries=result)
    finally:
        cursor.close()
        conn.close()


@app.route('/programs')
def programs():
    """list of programs"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute("""select p.name,count(distinct rp.regulon_id),count(distinct rpg.gene_id) from trans_programs p join regulon_programs rp on rp.program_id=p.id join regulon_program_genes rpg on rpg.program_id=p.id group by rp.program_id""")
        result = [{"program": program,
                   "num_regulons": num_regulons,
                   "num_genes": num_genes}
                  for program, num_regulons, num_genes in cursor.fetchall()]
        return jsonify(entries=result)
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


@app.route('/drug/<drugname>')
def drug(drugname):
    """information for the specified drug"""
    conn = dbconn()
    cursor = conn.cursor()
    try:
        cursor.execute("""
        select d.approved_symbol,d.approved,d.max_trial_phase,d.max_phase_gbm, dt.name as drug_type,at.name as action_type,ttm.name as target_type_model,moa.name as mechanism_of_action from drugs d join drug_types dt on d.drug_type_id=dt.id join action_types at on d.action_type_id=at.id join target_type_models ttm on d.target_type_model_id=ttm.id join mechanisms_of_action moa on d.mechanism_of_action_id=moa.id where d.name=%s""", [drugname])

        drugs = []
        for (approved_symbol, approved, max_trial_phase, max_phase_gbm,
             drug_type, action_type, target_type_model, mechanism_of_action
        ) in cursor.fetchall():
            drugs.append({
                "name": drugname,
                "approved": approved != 0,
                "max_trial_phase": max_trial_phase if max_trial_phase is not None else 0,
                "max_phase_gbm": max_phase_gbm if max_phase_gbm is not None else 0,
                "drug_type": drug_type,
                "action_type": action_type,
                "target_type_model": target_type_model,
                "mechanism_of_action": mechanism_of_action
            })
        return jsonify(drug=drugs[0])
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
        cursor.execute('select entrez_id,ensembl_id,preferred,uniprot_id,uniprot_function,ens_description from genes where entrez_id=%s or ensembl_id=%s or preferred=%s',
                       [gene, gene, gene])
        results = [row for row in cursor.fetchall()]
        if len(results) > 0:
            # not that "description" and "function" are currently swapped in the database
            entrez_id, ensembl_id, preferred, uniprot_id, description, function = results[0]
            return jsonify(entrez_id=entrez_id, ensembl_id=ensembl_id,
                           preferred=preferred, description=description,
                           uniprot_id=uniprot_id, function=function)
        else:
            return jsonify(entrez_id="NA", ensembl_id="NA", preferred="NA", description='NA',
                           uniprot_id='NA', function='NA')
    finally:
        cursor.close()
        conn.close()


if __name__ == '__main__':
    app.debug = True
    app.secret_key = 'trstrestnorgp654g'
    app.run(host='0.0.0.0', debug=True)
