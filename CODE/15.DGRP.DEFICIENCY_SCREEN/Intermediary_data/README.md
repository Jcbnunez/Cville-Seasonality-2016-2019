#README- description of intermediary files for Deficiency Analysis

1 . week1def.data, week2def.data,week3def.data : these are the accumulated behavior data from each of the three blocks of the assay.
Created using the corresponding "Data.collection" script, these files convert the trikinetic files from that block's measurements into a data form
that lists the activity for each fly, for each timepoint of collection.
2. datasmoothed90window30step.RDS: this file binds the three files above and aggregates activity over a sliding window.
Created using Data.aggregation.R, this file goes fly-by-fly, aggregating activity accross 90 second windows with 30 second steps
3. startle.data.RDS, this file gives the startle-duration phenototype and basal activity of each fly
Created using Phenotype.analysis.R.
4. emmeans.data.RDS, this file gives the activity-decay phenotype of each fly
Created using Phenotype.analysis.R.
5. modeling.data.RDS, this file gives the fly activity over time, to be used in creating complementation model tests
Created using Phenotype.analysis.R.