import os
from collections import defaultdict

def read_cancer_mutation(datadir, hr, hr_result, all_cancers, all_mutations):
    """Return a PMID dataset that is indexed for both cancers and mutations"""
    result = {}
    hr_result['cancer_mutation'] = result
    result['cancer'] = defaultdict(set)
    result['mutation'] = defaultdict(set)

    with open(os.path.join(datadir, 'Mut_Cancer_HR_%d.txt' % hr)) as infile:
        header = infile.readline().strip().split('\t')
        cancer_col = header.index('Disease')
        mutation_col = header.index('Mutation')
        pmids_col = header.index('PMIDs')
        for line in infile:
            row = line.strip().split('\t')
            if len(row) > 1:
                # synonyms for cancer
                cancers = [c for c in row[cancer_col].strip('"').split(' or ')
                           if len(c.strip()) > 0]
                mutation = row[mutation_col].strip()
                all_cancers.update(cancers)
                all_mutations.add(mutation)
                pmids = [s.strip() for s in row[pmids_col].split('|') if len(s.strip()) > 0]
                for c in cancers:
                    result['cancer'][c].update(pmids)
                for pmid in pmids:
                    hr_result['pmid_cancer'][pmid].update(cancers)
                    hr_result['pmid_mutation'][pmid].add(mutation)
                result['mutation'][mutation].update(pmids)


def read_mm_mutation(datadir, hr, hr_result, all_diseases, all_mutations):
    """Also known as disease-mutation"""
    result = {}
    hr_result['mm_mutation'] = result
    result['disease'] = defaultdict(set)
    result['mutation'] = defaultdict(set)

    with open(os.path.join(datadir, 'Mut_MM_HR_%d.txt' % hr)) as infile:
        header = infile.readline().strip().split('\t')
        disease_col = header.index('Disease')
        pmids_col = header.index('PMIDs')
        try:
            mutation_col = header.index('Mutation')
        except:
            print("WARNING: Dataset 'Mut_MM_HR_%d.txt' has no 'Mutation' column" % hr)
            mutation_col = header.index('Mutations')

        for line in infile:
            row = line.split('\t')
            diseases = [d for d in row[disease_col].strip('"').split(' or ') if len(d.strip()) > 0]
            all_diseases.update(diseases)
            mutations = [m.strip() for m in row[mutation_col].strip().split(' or ')
                         if len(m.strip()) > 0]
            if len(mutations) > 0:
                all_mutations.update(mutations)

            pmids = [s.strip() for s in row[pmids_col].split('|') if len(s.strip()) > 0]
            for d in diseases:
                result['disease'][d].update(pmids)
            for mutation in mutations:
                result['mutation'][mutation].update(pmids)

            for pmid in pmids:
                hr_result['pmid_disease'][pmid].update(diseases)
                hr_result['pmid_mutation'][pmid].update(mutations)


def read_mm_regulator(datadir, hr, hr_result, all_diseases, all_regulators):
    """also known as disease regulator"""
    result = {}
    hr_result['mm_regulator'] = result
    result['disease'] = defaultdict(set)
    result['regulator'] = defaultdict(set)

    with open(os.path.join(datadir, 'Regulator_MM_HR_%d.txt' % hr)) as infile:
        header = infile.readline().strip().split('\t')
        disease_col = header.index('Disease')
        regulator_col = header.index('Regulator')
        pmids_col = header.index('PMIDs')
        for line in infile:
            row = line.split('\t')
            if len(row) > 1:
                diseases = [d for d in row[disease_col].strip('"').split(' or ')
                            if len(d.strip()) > 0]
                all_diseases.update(diseases)
                regulator = row[regulator_col].strip()
                if len(regulator.strip()) > 0:
                    all_regulators.add(regulator)
                regulator_ensg = row[2].strip()
                pmids = [s.strip() for s in row[pmids_col].split('|') if len(s.strip()) > 0]
                for d in diseases:
                    result['disease'][d].update(pmids)
                result['regulator'][regulator].update(pmids)

                for pmid in pmids:
                    hr_result['pmid_disease'][pmid].update(diseases)
                    hr_result['pmid_regulator'][pmid].add(regulator)


def read_mm_regulon(datadir, hr, hr_result, all_diseases, all_regulons):
    """also known as disease regulon"""
    result = {}
    hr_result['mm_regulon'] = result
    result['disease'] = defaultdict(set)
    result['regulon'] = defaultdict(set)

    with open(os.path.join(datadir, 'Regulon_MM_HR_%d.txt' % hr)) as infile:
        header = infile.readline().strip().split('\t')
        disease_col = header.index('Disease')
        regulon_col = header.index('Regulon')
        pmids_col = header.index('PMIDs')
        for line in infile:
            row = line.split('\t')
            if len(row) > 1:
                diseases = [d for d in row[disease_col].strip('"').split(' or ')
                            if len(d.strip()) > 0]
                all_diseases.update(diseases)
                regulon = row[regulon_col]
                if len(regulon.strip()) > 0:
                    all_regulons.add(regulon)
                regulon_ensg = row[2]
                pmids = [s.strip() for s in row[pmids_col].split('|') if len(s.strip()) > 0]
                for d in diseases:
                    result['disease'][d].update(pmids)
                result['regulon'][regulon].update(pmids)

                for pmid in pmids:
                    hr_result['pmid_disease'][pmid].update(diseases)
                    hr_result['pmid_regulon'][pmid].add(regulon)


def read_mutation_regulator(datadir, hr, hr_result, all_mutations, all_regulators):
    """The input data sets are inconsistent !!! The Regulator column can either
    be anywhere or not there at all !"""
    result = {}
    hr_result['mutation_regulator'] = result
    result['mutation'] = defaultdict(set)
    result['regulator'] = defaultdict(set)

    with open(os.path.join(datadir, 'Mut_Regulators_Raw_HR_%d.txt' % hr)) as infile:
        header = infile.readline().strip().split('\t')
        mutation_col = header.index('Mutation')
        pmids_col = header.index('PMIDs')
        try:
            regulator_col = header.index('Regulator')
        except:
            # there is a multivalued Regulators column
            print("WARNING: no 'Regulator' column found in 'Mut_Regulators_Raw_HR_%d.txt'" % hr)
            regulator_col = header.index('Regulators')
        for line in infile:
            row = line.strip().split('\t')
            if len(row) > 1:
                mutation = row[mutation_col]
                if len(mutation.strip()) > 0:
                    all_mutations.add(mutation)
                else:
                    continue
                regulators = [r.strip() for r in row[regulator_col].strip().split(' or ')
                              if len(r.strip()) > 0]
                if len(regulators) > 0:
                    all_regulators.update(regulators)
                else:
                    continue
                pmids = [s.strip() for s in row[pmids_col].split('|') if len(s.strip()) > 0]
                result['mutation'][mutation].update(pmids)
                for regulator in regulators:
                    result['regulator'][regulator].update(pmids)

                for pmid in pmids:
                    hr_result['pmid_regulator'][pmid].update(regulators)
                    hr_result['pmid_mutation'][pmid].add(mutation)


def read_mutation_drugs(datadir, hr, hr_result, all_mutations, all_drugs):
    result = {}
    hr_result['mutation_drugs'] = result
    result['mutation'] = defaultdict(set)
    result['drug'] = defaultdict(set)

    with open(os.path.join(datadir, 'Mut_Drugs_Raw_HR_%d.txt' % hr)) as infile:
        header = infile.readline().strip().split('\t')
        mutation_col = header.index('Mutation')
        drugs_col = header.index('Drugs')
        pmids_col = header.index('PMIDs')
        for line in infile:
            row = line.split('\t')
            if len(row) > 1:
                mutation = row[mutation_col]
                drugs = [d for d in row[drugs_col].strip('"').split(' or ') if len(d.strip()) > 0]
                all_drugs.update(drugs)
                pmids = [s.strip() for s in row[pmids_col].split('|') if len(s.strip()) > 0]
                if len(mutation.strip()) > 0:
                    all_mutations.add(mutation)
                    result['mutation'][mutation].update(pmids)
                for drug in drugs:
                    result['drug'][drug].update(pmids)

                for pmid in pmids:
                    hr_result['pmid_drug'][pmid].update(drugs)
                    hr_result['pmid_mutation'][pmid].add(mutation)


def read_all_datasets(datadir):
    print('Reading in all datasets')
    result = {}
    all_cancers = set()
    all_diseases = set()
    all_mutations = set()
    all_regulators = set()
    all_regulons = set()
    all_drugs = set()

    for hr in range(2, 7):
        hr_result = {}
        # document assoctiations are global within an HR value
        hr_result['pmid_cancer'] = defaultdict(set)
        hr_result['pmid_mutation'] = defaultdict(set)
        hr_result['pmid_disease'] = defaultdict(set)
        hr_result['pmid_regulator'] = defaultdict(set)
        hr_result['pmid_regulon'] = defaultdict(set)
        hr_result['pmid_drug'] = defaultdict(set)

        result[hr] = hr_result
        print('Hazard ratio %d...' % hr)
        read_cancer_mutation(datadir, hr, hr_result, all_cancers, all_mutations)
        read_mm_mutation(datadir, hr, hr_result, all_diseases, all_mutations)
        read_mm_regulator(datadir, hr, hr_result, all_diseases, all_regulators)
        read_mm_regulon(datadir, hr, hr_result, all_diseases, all_regulons)
        read_mutation_regulator(datadir, hr, hr_result, all_mutations, all_regulators)
        read_mutation_drugs(datadir, hr, hr_result, all_mutations, all_drugs)

    result['cancers'] = sorted(all_cancers)
    result['diseases'] = sorted(all_diseases)
    result['mutations'] = sorted(all_mutations)
    result['regulators'] = sorted(all_regulators)
    result['regulons'] = sorted(all_regulons)
    result['drugs'] = sorted(all_drugs)
    return result


"""
Queries
"""


def disease_mutation(cur, hr, disease, mutation):
    # - join cancer and disease pmids
    dm_query = "select distinct pmid from exc_disease_mutation dm join exc_diseases d on dm.disease_id=d.id join exc_mutations m on dm.mutation_id=m.id where hr=%s"
    dm_params = [hr]
    if disease != 'All myelomas':
        dm_query += ' and d.name=%s'
        dm_params.append(disease)
    if mutation != 'All':
        dm_query += ' and m.name=%s'
        dm_params.append(mutation)
    cur.execute(dm_query, dm_params)
    dm_result = [row[0] for row in cur.fetchall()]

    cm_query = "select distinct pmid from exc_cancer_mutation cm join exc_cancers c on cm.cancer_id=c.id join exc_mutations m on cm.mutation_id=m.id where hr=%s"
    cm_params = [hr]

    if disease != 'All Cancers':
        cm_query += ' and c.name=%s'
        cm_params.append(disease)
    if mutation != 'All':
        cm_query += ' and m.name=%s'
        cm_params.append(mutation)

    cur.execute(cm_query, cm_params)
    cm_result = [row[0] for row in cur.fetchall()]
    return dm_result + cm_result


def disease_regulator(cur, hr, disease, regulator):
    query = "select distinct pmid from exc_disease_regulator dr join exc_diseases d on dr.disease_id=d.id join exc_regulators r on dr.regulator_id=r.id where hr=%s"
    params = [hr]
    if disease != 'All myelomas':
        query += ' and d.name=%s'
        params.append(disease)
    if regulator != 'All':
        query += ' and r.name=%s'
        params.append(regulator)
    cur.execute(query, params)
    return [row[0] for row in cur.fetchall()]


def disease_regulon(cur, hr, disease, regulon):
    query = "select distinct pmid from exc_disease_regulon dr join exc_diseases d on dr.disease_id=d.id join exc_regulons r on dr.regulon_id=r.id where hr=%s"
    params = [hr]
    if disease != 'All myelomas':
        query += ' and d.name=%s'
        params.append(disease)
    if regulon != 'All':
        query += ' and r.name=%s'
        params.append(regulon)
    cur.execute(query, params)
    return [row[0] for row in cur.fetchall()]


def mutation_regulator(cur, hr, mutation, regulator):
    query = "select distinct pmid from exc_mutation_regulator mr join exc_mutations m on mr.mutation_id=m.id join exc_regulators r on mr.regulator_id=r.id where hr=%s"
    params = [hr]
    if mutation != 'All':
        query += ' and m.name=%s'
        params.append(mutation)
    if regulator != 'All':
        query += ' and r.name=%s'
        params.append(regulator)
    cur.execute(query, params)
    return [row[0] for row in cur.fetchall()]


def mutation_drug(cur, hr, mutation, drug):
    query = "select distinct pmid from exc_mutation_drug md join exc_mutations m on md.mutation_id=m.id join exc_drugs d on md.drug_id=d.id where hr=%s"
    params = [hr]
    if mutation != 'All':
        query += ' and m.name=%s'
        params.append(mutation)
    if drug != 'All':
        query += ' and d.name=%s'
        params.append(drug)
    cur.execute(query, params)
    return [row[0] for row in cur.fetchall()]


def search_disease_mutation(cur, hr, search_term):
    dm_query = "select distinct pmid from exc_disease_mutation dm join exc_diseases d on dm.disease_id=d.id join exc_mutations m on dm.mutation_id=m.id where hr=%s and (d.name=%s or m.name=%s)"
    cur.execute(dm_query, [hr, search_term, search_term])
    dmres = [row[0] for row in cur.fetchall()]

    cm_query = "select distinct pmid from exc_cancer_mutation cm join exc_cancers c on cm.cancer_id=c.id join exc_mutations m on cm.mutation_id=m.id where hr=%s and (c.name=%s or m.name=%s)"
    cur.execute(cm_query, [hr, search_term, search_term])
    cmres = [row[0] for row in cur.fetchall()]
    return dmres + cmres


def search_disease_regulator(cur, hr, search_term):
    cur.execute("select distinct pmid from exc_disease_regulator dr join exc_diseases d on dr.disease_id=d.id join exc_regulators r on dr.regulator_id=r.id where hr=%s and (d.name=%s or r.name=%s)", [hr, search_term, search_term])
    return [row[0] for row in cur.fetchall()]


def search_disease_regulon(cur, hr, search_term):
    cur.execute("select distinct pmid from exc_disease_regulon dr join exc_diseases d on dr.disease_id=d.id join exc_regulons r on dr.regulon_id=r.id where hr=%s and (d.name=%s and r.name=%s)", [hr, search_term, search_term])
    return [row[0] for row in cur.fetchall()]


def search_mutation_regulator(cur, hr, search_term):
    cur.execute("select distinct pmid from exc_mutation_regulator mr join exc_mutations m on mr.mutation_id=m.id join exc_regulators r on mr.regulator_id=r.id where hr=%s and (m.name=%s or r.name=%s)", [hr, search_term, search_term])
    return [row[0] for row in cur.fetchall()]


def search_mutation_drug(cur, hr, search_term):
    cur.execute("select distinct pmid from exc_mutation_drug md join exc_mutations m on md.mutation_id=m.id join exc_drugs d on md.drug_id=d.id where hr=%s and (m.name=%s or d.name=%s)", [hr, search_term, search_term])
    return [row[0] for row in cur.fetchall()]
