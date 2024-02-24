CREATE TABLE acc2taxid (
        accession       VARCHAR(25) NOT NULL PRIMARY KEY,
        version         INTEGER,
        taxon_id        INTEGER,
        gi              INTEGER
    );
CREATE INDEX acc_idx ON acc2taxid (accession);
.mode tabs
.import /dev/stdin acc2taxid
