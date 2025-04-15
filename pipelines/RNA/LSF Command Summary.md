# LSF Command Summary for Running RNA-Seq Pipeline

LSF (**Load Sharing Facility**) is a **job scheduling system** used on HPC clusters to submit, monitor, and manage jobs efficiently. Below is a summary on how to use **LSF commands** for job submission, monitoring, and resource management.

---

## **📌 1. Submitting Jobs with `bsub`**
To submit a job to the LSF scheduler, use `bsub`:

```
bsub < run_pipeline.sh
```

or submit a command directly:

```
bsub -J myjob -o myjob.out -e myjob.err "echo Hello World"
```

**Common bsub Options**


| **Option** | **Description** |
|----------|----------------|
| -J jobname | Assigns a job name. | 
| -o output.txt | Saves standard output to output.txt. |
| -e error.txt | Saves error messages to error.txt. | 
| -n 8 | Requests 8 CPU cores. | 
| -W 4:00 | Sets a time limit of 4 hours. | 
| -R "rusage[mem=16000]" | Requests 16GB RAM. | 
| -q normal | Submits job to the normal queue. | 
	
Example: Submit a Job with 8 Cores and 16GB RAM

```
bsub -J RNAseq -o RNAseq.out -e RNAseq.err -n 8 -R "rusage[mem=16000]" -W 8:00 -q normal "bash run_pipeline.sh"
```


## **📌 2. Checking Job Status with `bjobs`**	

After submitting a job, use bjobs to check its status:

```
bjobs
```

Interpreting bjobs Output

```
JOBID   USER    STAT  QUEUE   FROM_HOST  EXEC_HOST  JOB_NAME   SUBMIT_TIME
123456  user1   RUN   normal  node01     node07     RNAseq     Dec 12 10:00
```

| **Column** | **Description**                                 |
|----------|-------------------------------------------------|
| JOBID  | Unique job ID.                                  | 
| USER  | Who submitted the job.                          |
| STAT | Job status (PEND, RUN, EXIT).                   | 
| QUEUE| Queue the job is in.                            | 
| FROM_HOST | Node where job is submitted.                    | 
| EXEC_HOST | Node where job is running.                      | 
| JOB_NAME | The name of the job.                            | 
| SUBMIT_TIME | The time Submiting the job. | 

**Check a Specific Job**
```
bjobs -l 123456
```

## **📌 3. Canceling or Suspending Jobs**	

Cancel a Job
```
bkill 123456
```
Cancel All Your Jobs
```
bkill 0
```
Suspend a Running Job
```
bstop 123456
```
Resume a Suspended Job
```
bresume 123456
```

## **📌 4. Checking Cluster Resources**	

Find Available Queues
```
bqueues
```
Check Resources on a Node
```
bhosts
```
Get Details of a Queue
```
bqueues -l normal
```
Find Maximum Memory & CPUs Available
To determine how much memory and CPUs you can request:
```
bqueues -l normal | grep "MAX MEM"
bqueues -l normal | grep "MAX NUM PROCESSORS"
```

To check real-time available resources:
```
bhosts -l
```
## **📌 4. 5. Running RNA-Seq Pipeline with LSF**	

To submit the RNA-seq pipeline from BioPipelineCraft, create an LSF script:
```
#!/bin/bash
#BSUB -J RNA_pipeline
#BSUB -o RNA_pipeline.%J.out
#BSUB -e RNA_pipeline.%J.err
#BSUB -n 16
#BSUB -R "rusage[mem=64000]"
#BSUB -W 24:00
#BSUB -q normal

# Load necessary modules
module load STAR
module load samtools
module load subread

# Run RNA-seq pipeline
bash run_pipeline.sh
```

Submit the job:
```
bsub < run_pipeline.sh
```


## **Summary of LSF Commands**

| **Command**                                       | **Purpose**                          |
|--------------------------------------------------|--------------------------------------|
| `bsub < script.sh`                               | Submit a job.                        |
| `bjobs`                                          | Check job status.                    |
| `bjobs -l 123456`                                | Detailed job info.                    |
| `bpeek 123456`                                   | View job output in real-time.        |
| `bkill 123456`                                   | Kill a specific job.                 |
| `bkill 0`                                       | Kill all jobs for the user.          |
| `bqueues`                                       | List all queues.                     |
| `bhosts`                                        | Check available nodes.               |
| `bqueues -l normal`                             | Get details of `normal` queue.       |




	
	
	
	
	
