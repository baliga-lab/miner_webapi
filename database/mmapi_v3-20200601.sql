DROP TABLE IF EXISTS `bc_mutation_tf`;
CREATE TABLE `bc_mutation_tf` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `bicluster_id` int(11) NOT NULL,
  `mutation_id` int(11) NOT NULL,
  `tf_id` int(11) NOT NULL,
  `role` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `bcmuttf_idx` (`bicluster_id`,`mutation_id`,`tf_id`)
) ENGINE=InnoDB AUTO_INCREMENT=27175 DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS `bc_mutation_tf_roles`;
CREATE TABLE `bc_mutation_tf_roles` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(50) NOT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=3 DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS `bc_tf`;
CREATE TABLE `bc_tf` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `bicluster_id` int(11) NOT NULL,
  `tf_id` int(11) NOT NULL,
  `role` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `bctf_idx` (`bicluster_id`,`tf_id`)
) ENGINE=InnoDB AUTO_INCREMENT=4263 DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS `bc_tf_roles`;
CREATE TABLE `bc_tf_roles` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(50) NOT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=3 DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS `bicluster_boxplot_data`;
CREATE TABLE `bicluster_boxplot_data` (
  `bicluster_id` int(11) DEFAULT NULL,
  `patient_id` int(11) DEFAULT NULL,
  `median` float DEFAULT NULL,
  `min_value` float DEFAULT NULL,
  `max_value` float DEFAULT NULL,
  `lower_quartile` float DEFAULT NULL,
  `upper_quartile` float DEFAULT NULL,
  UNIQUE KEY `bcbd_idx` (`bicluster_id`,`patient_id`),
  KEY `bbp_bicluster_id` (`bicluster_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS `bicluster_genes`;
CREATE TABLE `bicluster_genes` (
  `bicluster_id` int(11) NOT NULL,
  `gene_id` int(11) NOT NULL,
  UNIQUE KEY `bcgenes_idx` (`bicluster_id`,`gene_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS `bicluster_patients`;
CREATE TABLE `bicluster_patients` (
  `bicluster_id` int(11) NOT NULL,
  `patient_id` int(11) NOT NULL,
  UNIQUE KEY `bcpats_idx` (`bicluster_id`,`patient_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS `biclusters`;
CREATE TABLE `biclusters` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(50) NOT NULL,
  `cox_hazard_ratio` float DEFAULT NULL,
  `trans_program` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=15718 DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS `drugs`;
CREATE TABLE `drugs` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(200) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=135 DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS `exc_cancer_mutation`;
CREATE TABLE `exc_cancer_mutation` (
  `hr` int(11) DEFAULT NULL,
  `cancer_id` int(11) DEFAULT NULL,
  `mutation_id` int(11) DEFAULT NULL,
  `pmid` int(11) DEFAULT NULL,
  UNIQUE KEY `exc_can_mut_idx` (`hr`,`cancer_id`,`mutation_id`,`pmid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS `exc_cancers`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `exc_cancers` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=5 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;


DROP TABLE IF EXISTS `exc_disease_mutation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `exc_disease_mutation` (
  `hr` int(11) DEFAULT NULL,
  `disease_id` int(11) DEFAULT NULL,
  `mutation_id` int(11) DEFAULT NULL,
  `pmid` int(11) DEFAULT NULL,
  UNIQUE KEY `exc_dis_mut_idx` (`hr`,`disease_id`,`mutation_id`,`pmid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;


DROP TABLE IF EXISTS `exc_disease_regulator`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `exc_disease_regulator` (
  `hr` int(11) DEFAULT NULL,
  `disease_id` int(11) DEFAULT NULL,
  `regulator_id` int(11) DEFAULT NULL,
  `pmid` int(11) DEFAULT NULL,
  UNIQUE KEY `exc_dis_regul_idx` (`hr`,`disease_id`,`regulator_id`,`pmid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;


--
-- Table structure for table `exc_disease_regulon`
--

DROP TABLE IF EXISTS `exc_disease_regulon`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `exc_disease_regulon` (
  `hr` int(11) DEFAULT NULL,
  `disease_id` int(11) DEFAULT NULL,
  `regulon_id` int(11) DEFAULT NULL,
  `pmid` int(11) DEFAULT NULL,
  UNIQUE KEY `exc_dis_regulon_idx` (`hr`,`disease_id`,`regulon_id`,`pmid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

DROP TABLE IF EXISTS `exc_diseases`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `exc_diseases` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=22 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;


DROP TABLE IF EXISTS `exc_drugs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `exc_drugs` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=107 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;



--
-- Table structure for table `exc_mutation_drug`
--

DROP TABLE IF EXISTS `exc_mutation_drug`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `exc_mutation_drug` (
  `hr` int(11) DEFAULT NULL,
  `mutation_id` int(11) DEFAULT NULL,
  `drug_id` int(11) DEFAULT NULL,
  `pmid` int(11) DEFAULT NULL,
  UNIQUE KEY `exc_mut_drug_idx` (`hr`,`mutation_id`,`drug_id`,`pmid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;


DROP TABLE IF EXISTS `exc_mutation_regulator`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `exc_mutation_regulator` (
  `hr` int(11) DEFAULT NULL,
  `mutation_id` int(11) DEFAULT NULL,
  `regulator_id` int(11) DEFAULT NULL,
  `pmid` int(11) DEFAULT NULL,
  UNIQUE KEY `exc_mut_regulator_idx` (`hr`,`mutation_id`,`regulator_id`,`pmid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;


DROP TABLE IF EXISTS `exc_mutations`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `exc_mutations` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=135 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

DROP TABLE IF EXISTS `exc_regulators`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `exc_regulators` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=183 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;


--
-- Table structure for table `exc_regulons`
--

DROP TABLE IF EXISTS `exc_regulons`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `exc_regulons` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=265 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `genes`
--

DROP TABLE IF EXISTS `genes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `genes` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `ensembl_id` varchar(80) DEFAULT NULL,
  `entrez_id` varchar(80) DEFAULT NULL,
  `preferred` varchar(80) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `gene_ens_idx` (`ensembl_id`)
) ENGINE=InnoDB AUTO_INCREMENT=15093 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;


DROP TABLE IF EXISTS `hallmarks`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `hallmarks` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(200) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=19 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;


DROP TABLE IF EXISTS `lin_hallmarks`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `lin_hallmarks` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(200) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=22 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;


DROP TABLE IF EXISTS `mechanisms_of_action`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `mechanisms_of_action` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(200) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=7 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;


DROP TABLE IF EXISTS `mutations`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `mutations` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(50) NOT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=249 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

DROP TABLE IF EXISTS `patient_tf`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `patient_tf` (
  `patient_id` int(11) NOT NULL,
  `tf_id` int(11) NOT NULL,
  `tf_activity` float DEFAULT NULL,
  UNIQUE KEY `pattf_idx` (`patient_id`,`tf_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;


DROP TABLE IF EXISTS `patients`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `patients` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(100) NOT NULL,
  `pfs_survival` int(11) DEFAULT NULL,
  `pfs_status` int(11) DEFAULT NULL,
  `age` float DEFAULT NULL,
  `sex` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;


DROP TABLE IF EXISTS `regulon_drugs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `regulon_drugs` (
  `regulon_id` int(11) DEFAULT NULL,
  `drug_id` int(11) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;



DROP TABLE IF EXISTS `regulon_hallmarks`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `regulon_hallmarks` (
  `regulon_id` int(11) DEFAULT NULL,
  `hallmark_id` int(11) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;


DROP TABLE IF EXISTS `regulon_lin_hallmarks`;
CREATE TABLE `regulon_lin_hallmarks` (
  `regulon_id` int(11) DEFAULT NULL,
  `hallmark_id` int(11) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


--
-- Table structure for table `regulon_mechanism_of_action`
--

DROP TABLE IF EXISTS `regulon_mechanism_of_action`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `regulon_mechanism_of_action` (
  `regulon_id` int(11) DEFAULT NULL,
  `mechanism_of_action_id` int(11) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;


DROP TABLE IF EXISTS `regulon_target_class`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `regulon_target_class` (
  `regulon_id` int(11) DEFAULT NULL,
  `target_class_id` int(11) DEFAULT NULL,
  `pval` float DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

DROP TABLE IF EXISTS `target_classes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `target_classes` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(200) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=16 DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS `tfs`;
CREATE TABLE `tfs` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(50) NOT NULL,
  `cox_hazard_ratio` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `tfs_name_idx` (`name`)
) ENGINE=InnoDB AUTO_INCREMENT=479 DEFAULT CHARSET=latin1;


