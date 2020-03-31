=======
History
=======

0.11.0 (2019-12-31)
--------------------

* Operational wy2020 version


0.12.0 (TBD)
------------------

* Correct 'any' option for days of week report diagnostic option
* Limit density and precip calculations when database overwrite is False and those figures are not being generated
* Check topo.nc and snow.nc dimensions prior to processing
* SI units option
* Set default run_name field in CoreConfig.ini to avoid sql error
* Make compare_runs across water years
* Catch multiple database entries for a single field with exception when writing to the database doesn't exit cleanly
* Make output properties file names in write_properties more explicit by including date
* Change snow line calculation to percentile of dem with minimum snow amount
* Remove spatialnc and IPW file processing
* Single point values processing and independent database tables for time-series, single pixel input/output data
* Add elevation to stn_validate() figure
