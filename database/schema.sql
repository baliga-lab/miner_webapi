drop table if exists mutations;
drop table if exists regulons;
drop table if exists genes;
drop table if exists trans_programs;
drop table if exists regulon_regulator_roles;
drop table if exists regulon_mutation_regulator_roles;
drop table if exists regulon_genes;
drop table if exists regulon_programs;
drop table if exists regulon_mutation_regulator;
drop table if exists regulon_regulator;
drop table if exists regulon_program_genes;
drop table if exists drug_types;
drop table if exists action_types;
drop table if exists mechanisms_of_action;
drop table if exists drug_targets;
drop table if exists drug_regulons;
drop table if exists drug_programs;
drop table if exists drugs;
drop table if exists target_type_models;
drop table if exists cm_flow_types;
drop table if exists cmf_pathways;
drop table if exists cm_flows;

/* Pathway mutations */
create table mutations (id integer primary key auto_increment, name varchar(100) not null);
create table regulons (
  id integer primary key auto_increment,
  name varchar(50) not null,
  cox_hazard_ratio float
);
create table genes (
  id integer primary key auto_increment,
  ensembl_id varchar(80),
  entrez_id varchar(80),
  preferred varchar(80),
  is_mutation integer default 0,
  is_regulator integer default 0
);

/* create table tfs (id integer primary key auto_increment, name varchar(50) not null, symbol varchar(50)); */
create table trans_programs (id integer primary key auto_increment, name varchar(10) not null);
create table regulon_regulator_roles (id integer primary key auto_increment, name varchar(30) not null);
create table regulon_mutation_regulator_roles (id integer primary key auto_increment, name varchar(50) not null);

create table regulon_genes (regulon_id integer not null references regulons, gene_id integer not null references genes);
create table regulon_programs (regulon_id integer not null references regulons, program_id integer not null references trans_programs);

create table regulon_mutation_regulator (
  id integer primary key auto_increment,
  regulon_id integer not null references regulons,
  mutation_id integer not null references mutations,
  regulator_id integer not null references genes,
  role_id integer not null references regulon_mutation_regulator_roles
);

create table regulon_regulator (
  id integer primary key auto_increment,
  regulon_id integer not null references regulons,
  regulator_id integer not null references genes,
  role_id integer not null references regulon_regulator_roles
);

create table regulon_program_genes (
  id integer primary key auto_increment,
  regulon_id integer not null references regulons,
  program_id integer not null references trans_programs,
  gene_id integer not null references genes,
  is_disease_relevant integer
);

/*
 * Drug model
 */
create table drug_types (id integer primary key auto_increment, name varchar(200));

create table action_types (id integer primary key auto_increment, name varchar(200));

create table target_type_models (id integer primary key auto_increment, name varchar(200));

create table mechanisms_of_action (id integer primary key auto_increment, name varchar(200));

create table drug_targets (
  drug_id integer references drugs,
  gene_id integer references genes
);

create table drug_regulons (
  drug_id integer references drugs,
  regulon_id integer references regulons
);

create table drug_programs (
  drug_id integer references drugs,
  program_id integer references trans_programs
);

create table drugs (
  id integer primary key auto_increment,
  name varchar(200),
  approved_symbol varchar(200),
  approved integer,
  max_trial_phase integer,
  max_phase_gbm integer,
  drug_type_id integer not null references drug_types,
  action_type_id integer not null references action_types,
  target_type_model_id integer not null references target_type_models,
  mechanism_of_action_id integer not null references mechanisms_of_action
);


/*
 * Pre-populate database
 */
/*
Causal Mechanistic flow tables
Involved are the Mutation with a Regulator and a Regulon and associated values
The CM flow topology should be redundant with the regulator regulon mutation
associations and this is mostly to have more data and more efficient ways to
access the data
*/
create table cm_flow_types (id integer primary key auto_increment, name varchar(200));
create table cmf_pathways (id integer primary key auto_increment, name varchar(200));

create table cm_flows (
  id integer primary key auto_increment,
  cmf_id varchar(100),
  cmf_name varchar(500),
  cmf_type_id integer references cm_flow_types,
  mutation_id integer references mutations, /* can be null */
  cmf_pathway_id integer references cmf_pathways,
  mutation_gene_id integer references genes, /* can be null */
  regulator_id integer references genes,
  regulon_id integer references regulons,
  regulon_mutation_regulator_role_id integer references regulon_mutation_regulator_roles,
  regulon_mutation_regulator_pvalue float,
  regulon_regulator_spearman_r float,
  regulon_regulator_spearman_pvalue float,
  regulon_t_statistic float,
  regulon_log10_p_stratification float,
  fraction_edges_correctly_aligned float,
  fraction_aligned_diffexp_edges float,
  num_downstream_regulons integer,
  num_diffexp_regulons integer
);
