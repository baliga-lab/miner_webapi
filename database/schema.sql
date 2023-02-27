drop table if exists mutations;
drop table if exists biclusters;
drop table if exists genes;
drop table if exists tfs;
drop table if exists trans_programs;
drop table if exists bc_tf_roles;
drop table if exists bc_tf_bc_roles;
drop table if exists bc_mutation_tf_roles;
drop table if exists bicluster_genes;
drop table if exists bicluster_programs;
drop table if exists bc_mutation_tf;
drop table if exists bc_tf;
drop table if exists bc_tf_genes;
drop table if exists bc_program_genes;
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
create table biclusters (
  id integer primary key auto_increment,
  name varchar(50) not null,
  cox_hazard_ratio float
);
create table genes (
       id integer primary key auto_increment,
       ensembl_id varchar(80) not null,
       entrez_id varchar(80),
       preferred varchar(80),
       is_mutation integer default 0,
       is_regulator integer default 0
);

create table tfs (id integer primary key auto_increment, name varchar(50) not null, symbol varchar(50));
create table trans_programs (id integer primary key auto_increment, name varchar(10) not null);
create table bc_tf_roles (id integer primary key auto_increment, name varchar(30) not null);
create table bc_mutation_tf_roles (id integer primary key auto_increment, name varchar(50) not null);
create table bc_tf_bc_roles (id integer primary key auto_increment, name varchar(50) not null);

create table bicluster_genes (bicluster_id integer not null references biclusters, gene_id integer not null references genes);
create table bicluster_programs (bicluster_id integer not null references biclusters, program_id integer not null references trans_programs);

create table bc_mutation_tf (id integer primary key auto_increment, bicluster_id integer not null references biclusters, mutation_id integer not null references mutations, tf_id integer not null references tfs, role integer not null references bc_mutation_tf_roles);

create table bc_tf (
       id integer primary key auto_increment,
       bicluster_id integer not null references biclusters,
       tf_id integer not null references tfs,
       role integer not null references bc_tf_roles
);

create table bc_program_genes (
       id integer primary key auto_increment,
       bicluster_id integer not null references biclusters,
       program_id integer not null references tfs,
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
  regulon_id integer references biclusters
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
insert into bc_mutation_tf_roles (id, name) values (1, 'down-regulates');
insert into bc_mutation_tf_roles (id, name) values (2, 'up-regulates');

insert into bc_tf_roles (id, name) values (1, 'activates');
insert into bc_tf_roles (id, name) values (2, 'represses');

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
  tf_id integer references tfs,
  bc_id integer references biclusters,
  bc_mutation_tf_role_id integer references bc_mutation_tf_roles,
  bc_mutation_tf_pvalue float,
  bc_tf_spearman_r float,
  bc_tf_spearman_pvalue float,
  bc_t_statistic float,
  bc_log10_p_stratification float,
  fraction_edges_correctly_aligned integer,
  fraction_aligned_diffexp_edges float,
  num_downstream_regulons integer,
  num_diffexp_regulons integer
);
