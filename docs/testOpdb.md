## Run mcsActor writing to test opdb

The standard mcsActor can be run against a test opdb, where old PFSC data can be reprocessed and the measurements written to new tables, without overwriting any real opdb tables.

By default the mechanism needs to be run at Subaru using the real `mcsActor`. It is possible to run it otherwise, but that takes more fiddling.

### Create a test opdb

You first need to create a new test opdb. Well, actually, one of the admins needs to create one. Ask Yabe-san or Craig, and they will run `makeTestDb.sh $DBNAME`. We don't know if there are any interesting limits, but assume not. I suggest putting a username in the name; I have been using names like `cpl_mcs_01`. [To the admins, there are additional arguments to that command. In particular `-D` will `DROP` the database first and `-K` will kill the foreign server entirely. ]

Given that name and with `spt_operational_database` set up, `configTestDb.sh -s MCS $DBNAME` will configure that database such that all real opdb tables are readable, but the tables which are written by the `MCS` system are recreated (empty) in the test database. Specifically for now, `mcs_data`, `cobra_match`, `fiducial_fiber_match`.

The database runs on a server which is running on one of the Hilo servers (` pfsa-db01-gb.subaru.nao.ac.jp` on port 5434 (different from the usual 5432)).

For sanity, you will want something like the following in your `~/.pgpass`:

  `pfsa-db01-gb.subaru.nao.ac.jp:5434:*:pfs:PASSWORD`

### Set up the `mcsActor` and reprocess data

The `mcsActor` now tracks opdb connections a bit more carefully, and you can point it at the database of your choice with the `setDb` command. Given the test db as described above, use:

  `setDb hostname=pfsa-db01-gb.subaru.nao.ac.jp port=5434 db=cpl_mcs_01`

Sending an empty `setDb` command will reset the connection to the summit default.

After that, you can reprocess any old data on disk with, e.g. `expose object rerunFrameId=NNN doCentroid doFibreID`. That will read the PFSC for for that frame and the associated telescope information from the opdb, then re-centroid and re-id the spots, and write to the test tables. The PFSC file is not touched, nor are the real opdab tables. If rows already exist in the test db they are overwritten.

### Notes

 - For analysis and comparisons you should obviously make separate connections to the test and real `opdb`s. But for possible convenience, the real `opdb` tables are available in the test db as, say, `opdb_cobra_match`, etc.

 - we are not correctly versioning or conveniently linking the MCS boresight, so that will probably be wrong for rereductions. [Hmm, I thought we talked about adding a column to the `mcs_exposure` table.... -- CPL]

 - we can extend this to AGCC reductions easily enough, but have not yet.


 -------

 2024-04-05
