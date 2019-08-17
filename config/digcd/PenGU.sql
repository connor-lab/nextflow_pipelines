create table isolate(
	ID serial primary key not null,
	y_number text unique not null,
	episode_number text
    );

create table sequencing(
	ID SERIAL primary key not null,
	y_number text not null references isolate(y_number),
	sequencing_instrument text,
	sequencing_run text,
	sequencing_start_date date,
    sequencing_end_date date,
    calculated_magnitude int,
    trimmed_readlength decimal(5, 2),
    mean_insert_size int,
	qc_pass bool,
    created_at timestamp default current_timestamp
    );

create table mlst_scheme_metadata(
    ID SERIAL primary key not null,
    pubmlst_name text unique not null,
    pubmlst_url text unique not null,
    pubmlst_updated_at TIMESTAMP not null,
    locus_1 text,
    locus_2 text,
    locus_3 text,
    locus_4 text,
    locus_5 text,
    locus_6 text,
    locus_7 text
    );
    
create table mlst_sequence_types(
    ST text primary key unique not null,
    locus_1 text,
    locus_2 text,
    locus_3 text,
    locus_4 text,
    locus_5 text,
    locus_6 text,
    locus_7 text
    );

create table mlst(
    ID SERIAL primary key not null,
    y_number text not null references isolate(y_number),
    ST text references mlst_sequence_types(ST),
    locus_1 text,
    locus_2 text,
    locus_3 text,
    locus_4 text,
    locus_5 text,
    locus_6 text,
    locus_7 text,
    created_at timestamp default current_timestamp
    );

create table clustercode_snpaddress(
    ID SERIAL primary key not null,
    clustercode text unique not null,
    clustercode_frequency int,
    reference_name text not null,
    t250 int not null,
    t100 int not null,
    t50 int not null,
    t25 int not null,
    t10 int not null,
    t5 int not null,
    snpaddress_string text unique not null,
    created_at timestamp default current_timestamp,
    clustercode_updated timestamp not null
    );

create table clustercode(
    ID SERIAL primary key not null,
    y_number text not null references isolate(y_number),
    clustercode text not null references clustercode_snpaddress(clustercode),
    created_at timestamp default current_timestamp,
    clustercode_updated timestamp not null
    );

create table clustercode_history(
    ID SERIAL primary key not null,
    y_number text not null references isolate(y_number),
    old_clustercode text not null references clustercode_snpaddress(clustercode),
    new_clustercode text not null references clustercode_snpaddress(clustercode),
    created_at timestamp default current_timestamp
    );

create table reference_metadata(
    ID SERIAL primary key not null,
    reference_name text unique not null,
    snapperdb_db_name text unique not null,
    ncbi_accession text,
    genome_filename text unique not null,
    created_at timestamp default current_timestamp
    );

create table reference_distance(
    ID SERIAL primary key not null,
    y_number text not null references isolate(y_number),
    reference_name text not null references reference_metadata(reference_name),
    reference_mash_distance decimal(11, 10),
    reference_mash_p_value text not null,
    reference_common_kmers text not null,
    created_at timestamp default current_timestamp
    );

create table ribotype_metadata(
    ID SERIAL primary key not null,
    ribotype_reference_name text unique not null,
    ribotype text not null
    );

create table ribotype_distance(
    ID SERIAL primary key not null,
    y_number text not null references isolate(y_number),
    ribotype_reference_name text not null references ribotype_metadata(ribotype_reference_name),
    ribotype text not null,
    ribotype_mash_distance decimal(11, 10)
    );

create table toxinotype_db_tcda(
    ID SERIAL primary key not null,
    sequence text unique not null
    );

create table toxinotype_db_tcdb(
    ID SERIAL primary key not null,
    sequence text unique not null
    );

create table toxinotype_db_toxa(
    ID SERIAL primary key not null,
    sequence text unique not null
    );

create table toxinotype_db_toxb(
    ID SERIAL primary key not null,
    sequence text unique not null
    );

create table toxinotype(
    ID SERIAL primary key not null,
    y_number text not null references isolate(y_number),
    tox_a int references toxinotype_db_toxa(ID),
    tox_b int references toxinotype_db_toxb(ID),
    tcd_a int references toxinotype_db_tcda(ID),
    tcd_b int references toxinotype_db_tcdb(ID)
    );
