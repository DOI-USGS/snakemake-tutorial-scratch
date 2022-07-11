# Issue 4: Generalizing Rules with Wildcards

Pipeline steps: Compute the mean predicted temperature on each day of the year (at all depths) for each lake of interest

Concepts learned: wildcards, rule all, dry-run

## Update `calc_doy_means.py` script to use properties from Snakemake rule
Our current workflow is able to fetch and unzip raw lake temperature prediction data from ScienceBase. In the `out` file of our `1_fetch` step, we have a folder called `tmp` containing a set of csv files, each of which contains daily temperature predictions for a single lake (at multiple depths).

The daily data spans nearly 40 years, and we would like to process it by calculating the day-of-year mean temperature predictions (at all depths) for each lake. Let's add a new rule to do this.

For this tutorial, we will choose two lakes that we want to calculate the day-of-year mean temperatures for. The lakes that we want to process have the IDs `120020150` and `107072210`.

Let's process the data based on our current knowledge of Snakemake pipelining. We can write a rule that uses the `2_process/calc_doy_means.py` script. We will need to indicate the input and output files for the two lakes we have chosen in our rule as well:
```
rule calc_doy_means:
    input:
        "1_fetch/out/tmp/pgdl_nhdhr_120020150_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_107072210_temperatures.csv"
    output:
        "2_process/out/doy_120020150.csv",
        "2_process/out/doy_107072210.csv"
    script:
        "2_process/calc_doy_means.py"
```

This rule is perfectly acceptable and should work if you test it out by running `snakemake --cores 1 calc_doy_means` in your command line/Anaconda prompt.

Note: If you do test out the rule, be sure to delete the outputs before moving forward (`snakemake --cores 1 calc_doy_means --delete-all-output`). If you don't delete the outputs, the Snakemake pipeline will not execute next time because the outputs have already been generated!

Next, let's replace the hardcoded file names and lake IDs in the `2_process/calc_doy_means.py` script with our snakemake rule properties.

We have a multiple files as inputs and multiple files as outputs. We can treat our inputs and outputs as lists when we call them in Python.

To begin, we will be removing the hardcoded list of `lake_ids` in our script; however, we still need to indicate how many iterations our for loop needs to go through. We can choose the number of iterations by using the length of our input file list: `for i in range(len(snakemake.input))`.

Let's read in the `in_file` and `out_file` names for each for loop iteration. We can use the snakemake inputs and outputs property lists, and call the file name by index:
```
if __name__ == '__main__':
    for i in range(len(snakemake.input))
        out_file = snakemake.output[i]
        in_file = snakemake.input[i]
```

Lastly, let's extract the `lake_id` from our `out_file` name in each iteration of the for loop. Here's the full code block to update your `calc_doy_means.py` script with:
```
if __name__ == '__main__':
    for i in range(len(snakemake.input)):
        out_file = snakemake.output[i]
        in_file = snakemake.input[i]
        lake_id = os.path.splitext(out_file)[0].split('_')[-1]
        main(out_file, in_file, lake_id)
```

Try running this step again with `snakemake --cores 1 calc_doy_means` and make sure it works for you. If it does, delete the outputs before we move on with `snakemake --cores 1 calc_doy_means --delete-all-output`.

## Generalize your rule with wildcards
Snakemake allows us to generalize rules by using "wildcards".

Instead of listing an explicit filepath for the input and output file associated with each lake ID that we want to process, we can write just one input and one output with a wildcard in place of the varying `lake_id`. Wildcards are indicated by using `{}` around the wildcard name (`lake_id` in this case). Here's our updated rule that uses wildcards:
```
rule calc_doy_means:
    input:
        in_file = "1_fetch/out/tmp/pgdl_nhdhr_{lake_id}_temperatures.csv"
    output:
        out_file = "2_process/out/doy_{lake_id}.csv"
    script:
        "2_process/calc_doy_means.py"
```

## Resolve wildcards with `rule all`
We've set up our Snakefile with some wildcards, but we've lost the information about which `lake_ids` specifically we want to process. We would typically resolve our wildcards in the next pipeline rule as inputs, like this:
```
rule whatever_we_want_to_do_next:
    input:
        "2_process/out/doy_120020150.csv",
        "2_process/out/doy_107072210.csv"
```

This would allow Snakemake to intuit that the `lake_id` wildcard should have values `120020150` and `107072210`.

We aren't ready to build the next step of our pipeline yet, but we can still create another rule to resolve our wildcards. We could create a new rule with any name just like above, and run our pipeline by calling that named rule (in the example above, that would be `snakemake --cores 1 whatever_we_want_to_do_next`).

However, it is generally considered a good Snakemake practice to collect your final outputs in a `rule all`. The naming of this rule does not make any difference, but this is a standard convetion that is recommended by the source documentation.

The `rule all` is typically defined at the top of your Snakefile. The reason for placing it at the top of your Snakefile is that if no target rule is specified in the command line, Snakemake will use the first rule of the Snakefile as the target. This means that we could just run `snakemake --cores 1`, and Snakemake would use the first rule as your target. 

Let's add a `rule all` to the top of our Snakefile that will resolve our `lake_id` wildcard:
```
rule all:
    input:
        "2_process/out/doy_120020150.csv",
        "2_process/out/doy_107072210.csv"
```

## Test your Snakefile structure with a dry-run
Let's check to make sure we have built our wildcards correctly. So far, we've been testing our pipeline by running it and actually generating the outputs as we go. This can be a tedious way to test the pipeline because some steps may be slow to run.

We may want to test out our Snakefile to see if we have resolved all of our wildcards correctly without actually running the pipeline. We can do this by doing a `dry-run`. The `dry-run` command will allow you to see how the wildcard value is propagated to your input and output filenames in your command line. Try it out!
```
snakemake --dry-run 
```

## Update `calc_doy_means.py` script to use the new rule properties
Now that we know our wildcards can be correctly resolved from our Snakefile, we just have one more update to make.

Our `calc_doy_means.py` is currently set up to loop through our list of input files within the script. However, we have replaced our list of input files with wildcards. We want to update our `calc_doy_means.py` to expect a single input file. The input and output properties of our rule will be read into our python script with any wildcards already resolved, so we don't need to worry about formatting them with the `lake_id` in the script. We can also make use of the wildcards property of our rule to determine the `lake_id` instead of pullilng it out of a filename. Here's the code we want to update our `calc_doy_means.py` script with:
```
if __name__ == '__main__':
    out_file = snakemake.output['out_file']
    in_file = snakemake.input['in_file']
    lake_id = snakemake.wildcards['lake_id']
    main(out_file, in_file, lake_id)
```

We're all done - time to run our pipeline! Now that we have used wildcards in our Snakefile, that rule can actually be paralellized. Let's tell snakemake to use two cores so that each lake's data can be processed in parallel. Try it out:
```
snakemake --cores 2
```

We've successfully learned how to implement wildcards in our pipeline, and introducted the concepts of the `rule all` and how to test your pipeline with a `dry-run`.

In the next issue, we'll learn how to update your Snakefile with new wildcards, as well as how to combine these wildcards back together to produce a single file and plot for all of the lakes.