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
       preferred varchar(80)
);

create table tfs (id integer primary key auto_increment, name varchar(50) not null);
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

insert into bc_mutation_tf_roles (id, name) values (1, 'down-regulates');
insert into bc_mutation_tf_roles (id, name) values (2, 'up-regulates');

insert into bc_tf_roles (id, name) values (1, 'activates');
insert into bc_tf_roles (id, name) values (2, 'represses');
