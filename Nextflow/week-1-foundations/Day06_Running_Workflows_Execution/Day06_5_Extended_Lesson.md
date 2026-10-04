# Day 6: Running Workflows and Understanding Execution

**Learning Time**: 30 minutes  
**Prerequisites**: Days 1-5 completed, understanding of workflows  
**Goal**: Understand how Nextflow executes workflows and manages task execution

---

## 📖 Introduction (3 minutes)

Welcome to Day 6! Today you'll learn what happens **behind the scenes** when you run a Nextflow workflow. Understanding execution is crucial for debugging, optimizing, and trusting your pipelines.

**Why this matters**: Nextflow's execution model is what makes it powerful—isolated tasks, automatic parallelization, and resumability. Understanding how it works helps you debug problems, optimize performance, and appreciate why Nextflow handles 1,000 samples as easily as 1.

### What You'll Learn Today

- How to run Nextflow workflows (command syntax)
- What happens during execution
- The work directory structure and why it matters
- How resumability works (`-resume` is magic!)
- Reading and interpreting execution logs
- Monitoring workflow progress
- Understanding task isolation

### The Big Picture

When you write:
```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    results = PROCESS(samples)
}
```

Nextflow:
1. Creates an execution plan
2. Stages input files
3. Runs each task in isolation
4. Tracks what succeeded/failed
5. Collects outputs
6. Enables resume on failure

Today you'll understand each step!

---

## 🎯 Key Concepts with Examples (12 minutes)

### 1. Running a Workflow - Command Line Basics (3 minutes)

**Basic Execution**:
```bash
# Run a workflow
nextflow run workflow.nf

# Run with parameters
nextflow run workflow.nf --input data/ --output results/

# Resume a failed workflow
nextflow run workflow.nf -resume

# With execution reports
nextflow run workflow.nf -with-report report.html -with-trace trace.txt
```

**Command Structure**:
```bash
nextflow run <workflow.nf> [options] [--params]

Common options:
  -resume              Continue from last successful step
  -with-report FILE    Generate HTML execution report
  -with-trace FILE     Generate execution trace
  -with-timeline FILE  Generate execution timeline
  -with-dag FILE       Generate workflow DAG visualization
  -work-dir DIR        Specify work directory location
  -profile NAME        Use configuration profile
```

**Example with Multiple Options**:
```bash
nextflow run rnaseq_pipeline.nf \
  --input 'data/*.fastq' \
  --genome refs/hg38.fa \
  --outdir results/ \
  -resume \
  -with-report execution_report.html \
  -with-timeline timeline.html
```

**What Happens When You Run**:
```
$ nextflow run workflow.nf

N E X T F L O W  ~  version 23.10.0
Launching `workflow.nf` [peaceful_darwin] DSL2 - revision: 1a2b3c4d5e

executor >  local (12)
[a1/b2c3d4] process > FASTQC (1)        [100%] 4 of 4 ✔
[e5/f6g7h8] process > TRIM (1)          [100%] 4 of 4 ✔
[i9/j0k1l2] process > ALIGN (1)         [ 50%] 2 of 4

```

**Understanding the Output**:
- `[peaceful_darwin]`: Random run name for identification
- `executor > local (12)`: Running locally with 12 tasks total
- `[a1/b2c3d4]`: Work directory hash for this task
- `[100%]`: Progress percentage
- `4 of 4 ✔`: 4 completed successfully out of 4 total

### 2. The Work Directory - Where Magic Happens (3 minutes)

**Every task runs in its own isolated directory**:

```
work/
├── a1/
│   └── b2c3d4e5f6.../          # FASTQC task 1
│       ├── sample1.fastq       # Symlink to input
│       ├── sample1_fastqc.html # Output file
│       ├── sample1_fastqc.zip  # Output file
│       ├── .command.sh         # Script that was executed
│       ├── .command.run        # Wrapper script
│       ├── .command.out        # Standard output
│       ├── .command.err        # Standard error
│       ├── .command.log        # Execution log
│       ├── .command.begin      # Timestamp when started
│       └── .exitcode           # Exit status (0 = success)
│
├── e5/
│   └── f6g7h8i9j0.../          # TRIM task 1
│       ├── sample1.fastq
│       ├── sample1_trimmed.fastq
│       └── .command.*
│
└── [more task directories...]
```

**Key Files Explained**:

**`.command.sh`** - The actual script that ran:
```bash
#!/bin/bash
# Nextflow-generated script

fastqc sample1.fastq
```

**`.command.out`** - Standard output:
```
Started analysis of sample1.fastq
Approx 95% complete for sample1.fastq
Analysis complete for sample1.fastq
```

**`.command.err`** - Standard error (errors/warnings):
```
Warning: Duplicate read names found
```

**`.exitcode`** - Exit status:
```
0
```
(0 = success, non-zero = failure)

**Why Isolation Matters**:
1. **No conflicts**: Each task has its own directory
2. **Parallel safety**: Tasks can't interfere with each other
3. **Debugging**: Easy to inspect what happened
4. **Resumability**: Nextflow knows exactly what succeeded

**Inspecting a Task**:
```bash
# Navigate to a task directory
cd work/a1/b2c3d4e5f6.../

# See what script ran
cat .command.sh

# Check the output
cat .command.out

# Check for errors
cat .command.err

# See the exit code
cat .exitcode

# Run the exact command that failed
bash .command.run
```

### 3. Resumability - The Killer Feature (3 minutes)

**The Problem** (without Nextflow):
```python
# Python pipeline
for sample in samples:
    qc(sample)      # Takes 10 min
    align(sample)   # Takes 30 min
    variants(sample) # Takes 20 min

# If it crashes at sample 50/100:
# - Lost 50 hours of work!
# - Must restart from sample 1
# - No easy way to skip completed samples
```

**The Solution** (with Nextflow):
```bash
# First run - crashes at sample 50
nextflow run pipeline.nf
# [45 samples complete, 5 failed]

# Just add -resume!
nextflow run pipeline.nf -resume
# [skips 45 completed, retries 5 failed]
```

**How Resume Works**:

1. **Task Hashing**: Nextflow creates a unique hash for each task based on:
   - Process name
   - Input files (contents)
   - Script code
   - Container image
   - Parameters

2. **Cache Checking**: Before running, Nextflow checks:
   - "Have I seen this exact task before?"
   - "Did it complete successfully?"
   - "Are the outputs still there?"

3. **Smart Skipping**:
   - If hash matches and outputs exist: **SKIP** ✓
   - If hash doesn't match or outputs missing: **RUN** ⟳

**Example**:
```
$ nextflow run workflow.nf

[FASTQC] sample1 → Hash: a1b2c3... → Run → Success ✓
[FASTQC] sample2 → Hash: d4e5f6... → Run → Success ✓
[FASTQC] sample3 → Hash: g7h8i9... → Run → FAILED ✗

$ nextflow run workflow.nf -resume

[FASTQC] sample1 → Hash: a1b2c3... → Cached ✓ (skip)
[FASTQC] sample2 → Hash: d4e5f6... → Cached ✓ (skip)
[FASTQC] sample3 → Hash: g7h8i9... → Run again ⟳
```

**What Invalidates the Cache**:
- ✗ Changed script code in the process
- ✗ Different input file contents
- ✗ Different parameters
- ✗ Different container
- ✓ Same everything = cached!

**Resume Best Practices**:
```bash
# Always use -resume when rerunning
nextflow run pipeline.nf -resume

# Clean work directory occasionally (removes cache)
nextflow clean -f

# See what would be cached
nextflow log  # Shows previous runs
nextflow log <run_name> -f hash,status,name
```

### 4. Execution Logs - Your Debugging Friend (2 minutes)

**The `.nextflow.log` File**:

Every run creates a detailed log in `.nextflow.log`:

```bash
# View the log
less .nextflow.log

# Search for errors
grep ERROR .nextflow.log

# See what tasks failed
grep "Task failed" .nextflow.log
```

**Example Log Content**:
```
Jan-23 10:30:15.123 [main] INFO  nextflow.cli.Launcher - Launching pipeline.nf
Jan-23 10:30:15.456 [main] INFO  nextflow.Session - Starting session...
Jan-23 10:30:16.789 [Task monitor] DEBUG nextflow.processor.TaskRun - [a1/b2c3d4] Staging input file: sample1.fastq
Jan-23 10:30:17.012 [Task monitor] INFO  nextflow.processor.TaskRun - [a1/b2c3d4] Submitted task > FASTQC (1)
Jan-23 10:30:45.234 [Task monitor] DEBUG nextflow.processor.TaskRun - [a1/b2c3d4] Task completed > FASTQC (1)
Jan-23 10:30:45.345 [Task monitor] ERROR nextflow.processor.TaskRun - [e5/f6g7h8] Task failed > TRIM (1) -- Error: command failed with exit status: 1
```

**Understanding Log Levels**:
- `INFO`: Normal execution information
- `DEBUG`: Detailed technical information
- `WARN`: Warnings (not errors)
- `ERROR`: Something went wrong

**Finding Failed Tasks**:
```bash
# Get the work directory of failed tasks
grep "Task failed" .nextflow.log

# Example output:
# [e5/f6g7h8] Task failed > TRIM (1)

# Navigate to that directory
cd work/e5/f6g7h8*/

# Investigate
cat .command.err  # See the error
cat .command.sh   # See what ran
bash .command.run # Reproduce the error
```

### 5. Execution Reports - Visual Insights (1 minute)

**Generate Reports**:
```bash
nextflow run workflow.nf -with-report report.html
```

**What You Get**:

**HTML Report** (`report.html`):
- Summary statistics
- Resource usage (CPU, memory, time)
- Task duration distribution
- Success/failure breakdown

**Sample Report Content**:
```
Total Tasks: 156
Completed: 153
Failed: 3
Cached: 0

CPU Hours: 42.5
Peak Memory: 64 GB
Total Time: 2h 15m

Slowest Tasks:
1. ALIGN (sample_042): 45m
2. ALIGN (sample_087): 43m
3. VARIANTS (sample_023): 38m
```

**Execution Timeline** (`-with-timeline timeline.html`):
- Visual timeline of when tasks ran
- Shows parallelization
- Identifies bottlenecks

**Trace File** (`-with-trace trace.txt`):
- Detailed CSV of all tasks
- Perfect for analysis in Python/R
- Columns: task_id, status, exit, duration, cpu, memory, etc.

---

## 🔗 Python Connection: Script Execution vs Workflow Execution (3 minutes)

### Python Script Execution

```python
# run_pipeline.py
import subprocess
import sys

def main():
    samples = ["s1", "s2", "s3"]
    
    for sample in samples:
        try:
            # Run process
            result = subprocess.run(
                ["fastqc", f"{sample}.fastq"],
                check=True,
                capture_output=True
            )
            print(f"✓ Completed {sample}")
        except subprocess.CalledProcessError as e:
            print(f"✗ Failed {sample}: {e}")
            # What now? Continue? Stop? Retry?
            sys.exit(1)  # Usually just crash

if __name__ == "__main__":
    main()
```

**Execution**:
```bash
python run_pipeline.py
# ✓ Completed s1
# ✓ Completed s2
# ✗ Failed s3: Command failed with exit code 1

# To resume? Must manually track what succeeded!
```

**Issues**:
- No automatic tracking of success/failure
- No resume capability
- Manual error handling needed
- Hard to debug (where did output go?)
- No isolation (all in same directory)

### Nextflow Workflow Execution

```groovy
// workflow.nf
process FASTQC {
    input: val(sample)
    output: path("*.html")
    script: "fastqc ${sample}.fastq"
}

workflow {
    samples = Channel.from("s1", "s2", "s3")
    FASTQC(samples)
}
```

**Execution**:
```bash
nextflow run workflow.nf
# [a1/b2c3] FASTQC (s1) ✔
# [d4/e5f6] FASTQC (s2) ✔
# [g7/h8i9] FASTQC (s3) ✗

# Resume is built-in!
nextflow run workflow.nf -resume
# [a1/b2c3] FASTQC (s1) ✔ [cached]
# [d4/e5f6] FASTQC (s2) ✔ [cached]
# [g7/h8i9] FASTQC (s3) ⟳ [retrying]
```

**Benefits**:
- Automatic success/failure tracking
- Built-in resume
- Automatic error isolation
- Easy debugging (each task in own directory)
- Complete execution logs

### Comparison Table

| Aspect | Python Script | Nextflow Workflow |
|--------|--------------|-------------------|
| **Execution tracking** | Manual | Automatic |
| **Resume capability** | Manual checkpointing | `nextflow run -resume` |
| **Task isolation** | Same directory | Separate work dirs |
| **Debugging** | Check stdout/stderr | Inspect `.command.*` files |
| **Logs** | Manual logging | `.nextflow.log` |
| **Reports** | Manual | `-with-report` |
| **Parallelization** | Manual threading | Automatic |
| **Error recovery** | Manual retry logic | Automatic retry options |

---

## 💻 Hands-On Exercises (10 minutes)

### Exercise 1: Understanding Work Directory Structure (3 minutes)

**Scenario**: You ran this workflow and one task failed:

```groovy
process ANALYZE {
    input:
    val sample_id
    
    output:
    path "${sample_id}_results.txt"
    
    script:
    """
    analyze_tool ${sample_id} > ${sample_id}_results.txt
    """
}

workflow {
    Channel.from("sample1", "sample2", "sample3")
        | ANALYZE
}
```

**Output**:
```
[a1/b2c3d4] ANALYZE (sample1) ✔
[e5/f6g7h8] ANALYZE (sample2) ✗ Failed
[i9/j0k1l2] ANALYZE (sample3) ✔
```

**Questions**:
1. Where would you find the error message for sample2?
2. What files exist in the `work/e5/f6g7h8.../` directory?
3. How would you re-run just the failed command manually?

<details>
<summary>Click for answers</summary>

**1. Where to find the error**:
```bash
# Navigate to the failed task's work directory
cd work/e5/f6g7h8*/

# Check standard error
cat .command.err

# Or check the combined log
cat .command.log
```

**2. Files in the work directory**:
```
work/e5/f6g7h8.../
├── .command.sh         # The script: analyze_tool sample2 > sample2_results.txt
├── .command.run        # Wrapper script
├── .command.out        # Standard output (probably empty if it failed)
├── .command.err        # Standard error (contains the error message)
├── .command.log        # Combined log
├── .command.begin      # Timestamp when started
├── .exitcode           # Exit status (non-zero, indicating failure)
└── sample2_results.txt # May or may not exist, depending on failure type
```

**3. Re-run the failed command**:
```bash
# Navigate to the directory
cd work/e5/f6g7h8*/

# Option 1: Run the wrapper script (exactly as Nextflow did)
bash .command.run

# Option 2: Run the command directly
bash .command.sh

# Option 3: Run the actual command
analyze_tool sample2 > sample2_results.txt
```

This lets you:
- Test fixes before modifying the workflow
- Debug with additional flags
- See output in real-time
- Experiment with different parameters

**After fixing**: Update your workflow and run with `-resume` to skip successful tasks.
</details>

### Exercise 2: Resume Behavior Prediction (4 minutes)

**Initial Run**:
```bash
$ nextflow run pipeline.nf
[FASTQC] sample1 → Success ✔
[FASTQC] sample2 → Success ✔
[TRIM]   sample1 → Success ✔
[TRIM]   sample2 → Failed ✗
```

**For each scenario, predict what happens with `-resume`**:

**Scenario A**: You fix the TRIM process code and run:
```bash
nextflow run pipeline.nf -resume
```

**Scenario B**: You add a new parameter but don't change code:
```bash
nextflow run pipeline.nf -resume --quality_threshold 30
```

**Scenario C**: You replace `sample1.fastq` with a new version:
```bash
nextflow run pipeline.nf -resume
```

**Scenario D**: You just run `-resume` without any changes:
```bash
nextflow run pipeline.nf -resume
```

<details>
<summary>Click for answers</summary>

**Scenario A: Fixed TRIM code**
```
[FASTQC] sample1 → Cached ✔ (code didn't change)
[FASTQC] sample2 → Cached ✔
[TRIM]   sample1 → RE-RUN ⟳ (code changed, must re-run)
[TRIM]   sample2 → RE-RUN ⟳ (code changed, was failed anyway)
```
**Why**: Changed process code invalidates cache for that process

**Scenario B: Added parameter**
```
[FASTQC] sample1 → RE-RUN ⟳ (parameter affects hash)
[FASTQC] sample2 → RE-RUN ⟳
[TRIM]   sample1 → RE-RUN ⟳
[TRIM]   sample2 → RE-RUN ⟳
```
**Why**: New parameter changes task hash, even if process doesn't use it. Everything re-runs.

**Scenario C: Changed input file**
```
[FASTQC] sample1 → RE-RUN ⟳ (file content changed)
[FASTQC] sample2 → Cached ✔ (didn't change)
[TRIM]   sample1 → RE-RUN ⟳ (upstream changed)
[TRIM]   sample2 → RE-RUN ⟳ (was failed, retry)
```
**Why**: Input file hash changed, invalidating that sample's entire pipeline

**Scenario D: No changes**
```
[FASTQC] sample1 → Cached ✔
[FASTQC] sample2 → Cached ✔
[TRIM]   sample1 → Cached ✔
[TRIM]   sample2 → RE-RUN ⟳ (was failed, retry)
```
**Why**: Failed tasks always retry, successful tasks cached

**The Pattern**:
- Failed tasks → Always retry
- Successful tasks → Cache if nothing changed
- Changes → Invalidate cache for affected tasks
</details>

### Exercise 3: Reading Execution Traces (3 minutes)

**Given this trace file excerpt** (simplified):

```csv
task_id,hash,name,status,exit,duration,realtime,cpu,%cpu,memory
1,a1/b2c3d4,FASTQC (sample1),COMPLETED,0,5m 30s,5m 15s,4.5h,85.7%,2.1 GB
2,e5/f6g7h8,FASTQC (sample2),COMPLETED,0,5m 45s,5m 30s,4.6h,83.3%,2.2 GB
3,i9/j0k1l2,TRIM (sample1),COMPLETED,0,3m 20s,3m 10s,2.8h,88.4%,1.5 GB
4,m3/n4o5p6,TRIM (sample2),FAILED,1,1m 5s,1m 0s,0.9h,90.0%,1.1 GB
5,q7/r8s9t0,ALIGN (sample1),COMPLETED,0,25m 15s,24m 50s,23h,92.5%,8.3 GB
```

**Questions**:
1. Which task used the most memory?
2. Which task was most CPU-efficient?
3. Why did task 4 fail (based on available info)?
4. Which task took the longest wall-clock time?
5. If you had 4 CPUs, how long would tasks 1-3 take?

<details>
<summary>Click for answers</summary>

**1. Most memory**: Task 5 (ALIGN sample1) used 8.3 GB

**2. Most CPU-efficient**: Task 4 (TRIM sample2) at 90.0% CPU usage
- Despite failing, it used CPU efficiently
- Task 5 at 92.5% is slightly higher

**3. Why task 4 failed**:
- Exit code = 1 (non-zero = failure)
- Very short runtime (1m 5s) suggests it failed early
- Possible reasons: input validation failed, missing dependencies, corrupted input
- To investigate: `cd work/m3/n4o5p6*/ && cat .command.err`

**4. Longest wall-clock time**: Task 5 (ALIGN) at 25m 15s
- Even though FASTQC tasks had more CPU hours, ALIGN took longest real time

**5. With 4 CPUs, tasks 1-3**:
- Tasks 1 and 2 run in parallel: max(5m 30s, 5m 45s) = 5m 45s
- Task 3 starts after: 5m 45s + 3m 20s = 9m 5s total
- **Sequential would be**: 5m 30s + 5m 45s + 3m 20s = 14m 35s
- **Parallel saves**: 14m 35s - 9m 5s = 5m 30s (38% faster)

**Key insights**:
- CPU hours ≠ real time (parallel execution)
- Exit code tells you if task succeeded
- Memory and CPU usage help optimize resource allocation
- Short failed tasks often mean early validation errors
</details>

---

## 🤔 Reflection Activity (4 minutes)

### Question 1: The Value of Resume

**Scenario**: You're running a variant calling pipeline on 500 whole-genome samples. Each sample takes about 6 hours to process through the entire pipeline (QC → Align → Variants).

The pipeline crashes after 2 days when 350 samples are complete.

**Questions**:
1. Without `-resume`, how much time would you lose?
2. With `-resume`, what happens?
3. What's the total time saved?
4. Why is this critical for large-scale bioinformatics?

<details>
<summary>Click for thoughts</summary>

**1. Time lost without resume**:
- 350 samples × 6 hours = 2,100 hours of computation
- All wasted! Must start from zero
- Another 2 days (48 hours) to get back to where you were
- Total: 4 days to complete what should have taken 2 days

**2. With resume**:
- Nextflow skips the 350 successful samples
- Only processes the remaining 150 samples
- Continues exactly where it left off
- No wasted computation

**3. Time saved**:
- Without resume: 4 days total
- With resume: 2 days (original) + time for 150 samples ≈ 2.9 days
- **Saved: ~1.1 days of compute time**
- More importantly: **2,100 CPU hours not wasted**

**4. Why critical**:
- Large datasets are common in genomics
- Failures happen (network issues, memory errors, cluster maintenance)
- Time is expensive (researcher time, compute costs)
- Reproducibility: can verify results without re-running everything
- Experimentation: can tweak downstream steps without re-running upstream

**Real-world impact**:
- Saves days/weeks of compute time
- Reduces cloud computing costs ($$$$)
- Enables iterative development
- Makes large-scale analysis practical
</details>

### Question 2: Debugging Workflow

You run a workflow and get this error:

```
[ERROR] Process `CALL_VARIANTS (sample_042)` terminated with an error exit status (137)
```

**Questions**:
1. What does exit code 137 typically mean?
2. Where would you look to investigate?
3. What would you check first?
4. How would you fix it and continue?

<details>
<summary>Click for answers</summary>

**1. Exit code 137 meaning**:
- Exit code 137 = 128 + 9
- Signal 9 = SIGKILL (killed by system)
- **Typical cause: Out of memory (OOM)**
- Process was killed because it used too much RAM

**2. Where to investigate**:
```bash
# Find the work directory from the log
grep "sample_042" .nextflow.log
# Output shows: [e5/f6g7h8] Task failed...

# Navigate there
cd work/e5/f6g7h8*/

# Check what happened
cat .command.err  # Might show "Killed" or memory error
cat .command.log  # Full log
cat .exitcode     # Confirms 137
```

**3. What to check first**:
1. **Memory usage**: Was the process allocated enough memory?
2. **Input size**: Is sample_042 unusually large?
3. **Tool behavior**: Does variant caller load entire BAM in memory?
4. **System resources**: Was cluster node full?

**4. How to fix and continue**:

**Option A**: Increase memory for the process
```groovy
process CALL_VARIANTS {
    memory '16 GB'  // Increase from default
    
    // or dynamic:
    memory { 8.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3
    
    // ... rest of process
}
```

**Option B**: Just for this sample (parameter)
```bash
nextflow run pipeline.nf -resume --max_memory 32GB
```

**Then resume**:
```bash
nextflow run pipeline.nf -resume
```

All successful samples are cached, only sample_042 (and any downstream) re-run with more memory.

**Prevention**: Always specify appropriate memory in process directives!
</details>

### Question 3: Work Directory Management

Your workflow has run successfully many times over several weeks. Your `work/` directory is now 500 GB!

**Questions**:
1. Is it safe to delete the work directory?
2. What happens to `-resume` if you delete it?
3. How can you clean it safely?
4. What should you keep, what can you delete?

<details>
<summary>Click for answers</summary>

**1. Safe to delete?**
- **Yes**, if you've saved important results to `publishDir`
- **No**, if you might need to resume
- The work directory is **temporary** execution space

**2. What happens to resume?**
- Resume will **not work** - no cached results
- Workflow will run from scratch
- All tasks must re-execute
- It's like running for the first time

**3. How to clean safely**:

**Option 1**: Delete old runs (keep recent)
```bash
# List all runs
nextflow log

# Clean specific old runs
nextflow clean -n  # Dry run (see what would be deleted)
nextflow clean -f  # Force delete

# Or by date
find work/ -type d -mtime +30 -exec rm -rf {} +  # Older than 30 days
```

**Option 2**: Clean only failed tasks
```bash
# Nextflow can clean failed task directories
nextflow clean -f -k
```

**Option 3**: Complete cleanup
```bash
# Remove entire work directory
rm -rf work/

# Next run will be fresh
nextflow run pipeline.nf  # No -resume possible
```

**4. What to keep/delete**:

**Keep**:
- ✓ Published results (`publishDir` outputs)
- ✓ Important reports/logs
- ✓ Recent runs (last 1-2 weeks) if still developing

**Can delete**:
- ✗ Old work directories (>30 days)
- ✗ Failed task directories (if not debugging)
- ✗ Runs from old versions of workflow
- ✗ Intermediate files (if finals are published)

**Best practice**:
```groovy
process IMPORTANT_STEP {
    publishDir "results/", mode: 'copy'  // Keep important outputs
    
    input:
    path input
    
    output:
    path "important_result.txt"
    
    // work/ directory can be cleaned
    // results/ directory is permanent
}
```

**Automation**:
```bash
# Add to cron job
0 0 * * 0 nextflow clean -f -before '1 month ago'
```
</details>

---

## 📝 Key Takeaways

Before moving to Day 7, ensure you understand:

✅ **How to run workflows** from command line with options  
✅ **The work directory structure** and why tasks are isolated  
✅ **Resume works by task hashing** and checking cache  
✅ **`.nextflow.log` contains** execution details for debugging  
✅ **Work directory files** (`.command.sh`, `.command.err`, etc.)  
✅ **Execution reports** provide insights into performance  
✅ **Task isolation enables** parallel, resumable execution  

### The Mental Model: Execution Flow

```
1. Parse workflow.nf
     ↓
2. Create execution plan (DAG)
     ↓
3. For each task:
   - Calculate hash
   - Check cache
   - If cached: skip
   - If not: run in work/xx/yyyy/
     ↓
4. Monitor execution
     ↓
5. Log everything
     ↓
6. Collect outputs
```

### Essential Commands Reference

```bash
# Running
nextflow run workflow.nf                    # Basic execution
nextflow run workflow.nf -resume            # Resume from failures
nextflow run workflow.nf --param value      # With parameters

# Reporting
nextflow run workflow.nf -with-report report.html
nextflow run workflow.nf -with-trace trace.txt
nextflow run workflow.nf -with-timeline timeline.html
nextflow run workflow.nf -with-dag dag.html

# Debugging
nextflow log                                # List runs
nextflow log <run_name> -f hash,status      # Task details
cd work/xx/yyyy*/                           # Navigate to task
cat .command.err                            # See error

# Cleanup
nextflow clean -n                           # Preview cleanup
nextflow clean -f                           # Delete work dirs
nextflow clean -f -k                        # Keep latest
```

---

## 🎯 Ready for Day 7?

Tomorrow is Day 7 - **Week 1 Review and Integration**! You'll consolidate everything you've learned this week by building a complete project that uses all concepts from Days 1-6.

### Week 1 Summary - What You've Learned

✅ **Day 1**: What Nextflow is and why it matters  
✅ **Day 2**: Groovy essentials for Nextflow  
✅ **Day 3**: Writing Nextflow processes  
✅ **Day 4**: Understanding channels  
✅ **Day 5**: Connecting processes into workflows  
✅ **Day 6**: Running and debugging workflows  

**Tomorrow**: Integrate everything into a complete quality control pipeline!

---

## ✅ Day 6 Completion Checklist

Before marking Day 6 complete, ensure you can:

- [ ] Understand `nextflow run` command syntax
- [ ] Explain the work directory structure
- [ ] Describe how resume works (task hashing)
- [ ] Know where to find error messages
- [ ] Understand what `.command.*` files contain
- [ ] Explain task isolation benefits
- [ ] Use execution reports for debugging

**Self-Test**: You have a failed workflow. You need to:
1. Find which task failed
2. See the error message
3. Re-run just that task manually
4. Fix it and resume

Can you describe the steps?

<details>
<summary>Check your answer</summary>

**Steps to debug and fix**:

1. **Find failed task**:
   ```bash
   # Check the console output or log
   grep "Task failed" .nextflow.log
   # Shows: [e5/f6g7h8] Task failed > PROCESS (sample)
   ```

2. **See error message**:
   ```bash
   cd work/e5/f6g7h8*/
   cat .command.err
   ```

3. **Re-run manually**:
   ```bash
   # Still in work/e5/f6g7h8*/
   bash .command.run
   # Or: bash .command.sh
   ```

4. **Fix and resume**:
   ```bash
   # Edit workflow.nf to fix the issue
   vim workflow.nf
   
   # Resume - successful tasks cached, failed ones retry
   nextflow run workflow.nf -resume
   ```

If you got this right, you can debug Nextflow workflows! 🎉
</details>

**Completed Day 6?** Update your `PROGRESS.md`! You understand Nextflow execution! 🚀

**Your progress**: 6/28 days (21.4%) complete

---

*Tomorrow: Day 7 - Week 1 Review and Integration Project*

**Almost done with Week 1! See you tomorrow! 🚀**