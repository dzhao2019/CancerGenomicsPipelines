# Day 1 Quick Reference Card

Print this or bookmark it for quick lookup while learning!

---

## 🎯 Day 1 Learning Objectives: Did You Achieve Them?

- [ ] Articulate the problem (why Python scripts become difficult)
- [ ] Describe what Nextflow is (orchestration tool)
- [ ] Identify scenarios where Nextflow provides value
- [ ] Understand when to use Python vs Nextflow
- [ ] Grasp the 3 core concepts: processes, channels, workflows
- [ ] Appreciate Nextflow's philosophy: reproducibility, portability, scalability

---

## 🔑 Core Concepts at a Glance

### Process (What Happens)

```groovy
process DO_SOMETHING {
    input:
        path input_file
    output:
        path "*.result"
    script:
        """
        some_command ${input_file} > output.result
        """
}
```

**Key points:**
- Isolated unit of work
- Can run anywhere (laptop, cluster, cloud)
- Can run many times in parallel
- Input/output explicitly declared

### Channel (How Data Flows)

```groovy
// Create a channel
samples = Channel.fromPath("data/*.fastq")

// Channels automatically parallelize
process_files(samples)  // Runs once for EACH file
```

**Key points:**
- Like a conveyor belt of data
- Multiple items = automatic parallelization
- Pass between processes
- Two types: queue (used once), value (reusable)

### Workflow (How It's Orchestrated)

```groovy
workflow {
    input = Channel.fromPath("data/*.fastq")
    out1 = PROCESS_1(input)
    out2 = PROCESS_2(out1)
}
```

**Key points:**
- Defines which processes run
- How data flows between them
- Reads declaratively (like a narrative)
- Nextflow handles scheduling and resources

---

## 🐍 Python vs Nextflow Quick Lookup

### Single Item Processing

**Python:**
```python
data = read_file("sample.fastq")
result = analyze(data)
```

**When:** Single file, data science work  
**Tools:** pandas, numpy, scipy, matplotlib

### Multiple Items Sequential

**Python:**
```python
for sample in samples:
    result = process(sample)
```

**When:** Few samples (<10), simple logic  
**Speed:** Slow (hours/days)

### Multiple Items Parallel (Python)

```python
from multiprocessing import Pool

with Pool(4) as p:
    results = p.map(process, samples)
```

**When:** Parallelization needed  
**Problem:** Complex code, no resumability

### Multiple Items Parallel (Nextflow)

```groovy
samples = Channel.fromPath("*.fastq")
process_files(samples)  // Automatic parallelization
```

**When:** Multiple samples, multiple tools, production pipelines  
**Speed:** Fast (minutes/hours)

---

## ⏱️ Time Comparison Examples

### Scenario: 1000 samples, 50 minutes per sample

| Approach | Time | Notes |
|----------|------|-------|
| Python `for` loop | 833 hours (35 days) | Sequential |
| Python `multiprocessing` (4 cores) | 208 hours (8.7 days) | Manual parallel |
| Nextflow (8 cores) | 106 hours (4.4 days) | Automatic parallel |
| Nextflow (64 cores cluster) | 13 hours | Scales automatically |

### Scenario: Pipeline fails at sample 500/1000

| Approach | Recovery | Total Time |
|----------|----------|-----------|
| Python `for` | Manual restart from 1 | +833 hours |
| Nextflow | `nextflow run -resume` | +4 hours (continue) |

---

## 📋 Python → Nextflow Decision Tree

```
┌─ Is this data processing / statistics?
│  └─ YES → Use Python + libraries
│
└─ Need to coordinate multiple tools?
   ├─ YES + many samples?
   │  ├─ YES → Use Nextflow (massive speedup)
   │  └─ NO → Could be either (Python simpler)
   │
   └─ NO → Use Python
```

**Simple rules:**

1. **One file + statistics/science** → Python
2. **Many files + many tools** → Nextflow
3. **Many files + Python processing** → Nextflow calling Python
4. **Unsure** → Ask "Will this run on 1000 samples?"

---

## 🚨 Common Misconceptions CORRECTED

| Wrong | Right |
|-------|-------|
| "Nextflow replaces Python" | Nextflow orchestrates; Python processes. Different tools for different jobs. |
| "Nextflow is a programming language" | Nextflow is a workflow orchestration system (uses Groovy, bash, Python, etc.) |
| "Channels are like Python lists" | Channels are data streams; lists are in-memory collections. Very different! |
| "Processes are like functions" | Processes are isolated units that can run anywhere; functions are local |
| "Use Nextflow for everything" | Use the right tool: Python for data science, Nextflow for orchestration |

---

## 🔄 The Nextflow Advantage: Resumability

### Without Nextflow (Python)
```python
# Pipeline fails at sample 500
for i, sample in enumerate(samples):  # 1000 samples
    try:
        process_sample(sample)  # FAILS at sample 500
    except:
        save_checkpoint()
        # Now you have to manually restart from 500
        # How do you know which succeeded? Is data consistent?
```

### With Nextflow
```
First run:
$ nextflow run pipeline.nf
# Fails at sample 500

Second run (resume):
$ nextflow run pipeline.nf -resume
# Continues from sample 501
# Automatic! No manual checkpointing!
```

**Time saved:** If each sample takes 50 min, first 499 samples × 50 = 24,950 minutes
- Python: Run all 50,000 minutes again
- Nextflow: Continue from 500 (saves 24,950 minutes!)

---

## 💾 File Processing Comparison

### Python Approach
```python
# Load everything into memory
import glob
files = glob.glob("data/*.fastq")  # All files into list
for f in files:
    data = load(f)  # Load one at a time
    result = process(data)
    # Risk: If memory-intensive, will crash with many files
```

**Problems:**
- Must know file list upfront
- Everything sequential by default
- Memory-intensive with large files

### Nextflow Approach
```groovy
// Stream files as they're processed
samples = Channel.fromPath("data/*.fastq")
process_files(samples)
// Nextflow handles scheduling and memory automatically
```

**Advantages:**
- Discovers files dynamically
- Parallelizes automatically
- Memory-efficient
- Works with 1 or 1,000,000 files

---

## 🛠️ Tool Coordination Comparison

### Python
```python
# You write all this manually:
import subprocess

result1 = subprocess.run(["tool1", "input.txt"])
result2 = subprocess.run(["tool2", result1.stdout])
result3 = subprocess.run(["tool3", result2.stdout])

# Problems:
# - Error handling?
# - What if tool2 fails?
# - What if you have 1000 inputs?
# - Parallelization?
```

### Nextflow
```groovy
process TOOL1 {
    script: "tool1 ${input}"
}
process TOOL2 {
    script: "tool2 ${input}"
}
process TOOL3 {
    script: "tool3 ${input}"
}

workflow {
    files = Channel.fromPath("*.txt")
    t1 = TOOL1(files)
    t2 = TOOL2(t1)
    t3 = TOOL3(t2)
}

// Automatic:
// - Parallelization
// - Error handling
// - Resumability
// - Logging
```

---

## 📚 Nextflow Concepts Reference

### Process Block Components

```groovy
process PROCESS_NAME {
    // 1. Directives (instructions to Nextflow)
    cpus 4
    memory '8GB'
    time '1h'
    publishDir "results"
    
    // 2. Input (what does it need)
    input:
        path reads
        val sample_id
    
    // 3. Output (what does it produce)
    output:
        path "*.bam"
    
    // 4. Script (the actual work)
    script:
        """
        bwa mem reference.fa ${reads} > output.bam
        """
}
```

### Common Directives

| Directive | Purpose | Example |
|-----------|---------|---------|
| `cpus` | CPU cores needed | `cpus 8` |
| `memory` | RAM needed | `memory '16GB'` |
| `time` | Maximum runtime | `time '2h'` |
| `container` | Docker/Singularity image | `container 'biocontainers/bwa:latest'` |
| `publishDir` | Where to save outputs | `publishDir "results"` |
| `label` | Group processes | `label 'short_running'` |

### Common Input Types

| Type | Use Case | Example |
|------|----------|---------|
| `path` | File input | `path reads` |
| `val` | Simple value | `val sample_id` |
| `tuple` | Multiple related items | `tuple val(id), path(file)` |
| `stdin` | Standard input | `stdin` |

### Common Output Types

| Type | Use Case | Example |
|------|----------|---------|
| `path` | Output files | `path "*.bam"` |
| `stdout` | Console output | `stdout` |
| `val` | Simple value | `val "done"` |

---

## 🏃 Quick Start Template

```groovy
params {
    input_dir = "data/"
    output_dir = "results/"
}

process MY_PROCESS {
    input:
        path input_file
    output:
        path "*.output"
    script:
        """
        my_tool ${input_file} > output.txt
        """
}

workflow {
    // Find input files
    files = Channel.fromPath("${params.input_dir}/*.fastq")
    
    // Run process
    results = MY_PROCESS(files)
    
    // Optional: view results
    results.view()
}
```

---

## 🎓 Day 1 Completion Checklist

**Understanding Concepts:**
- [ ] Nextflow is orchestration (not programming)
- [ ] Python is data processing
- [ ] Together they're powerful
- [ ] Processes are isolated units
- [ ] Channels enable parallelization

**Appreciating Benefits:**
- [ ] Understand resumability value (saved ~25,000 minutes!)
- [ ] See automatic parallelization benefit
- [ ] Recognize portability advantage
- [ ] Know when Nextflow is overkill (single file)

**Ready for Day 2:**
- [ ] Know I need to learn Groovy (not deeply, just enough)
- [ ] Excited to write first process
- [ ] Understand this investment will pay off

---

## 💡 Key Quote to Remember

> "Python processes the data. Nextflow orchestrates the tools. Together, they're unstoppable."

---

## 🔗 Fast Reference Links

- **Official Docs:** https://www.nextflow.io/docs/latest/
- **Nextflow Patterns:** https://nextflow-io.github.io/patterns/
- **nf-core (Community):** https://nf-co.re/
- **Slack Community:** https://nf-core.slack.com

---

## 📝 Notes for Tomorrow (Day 2)

Remember these when learning Groovy:

1. **String interpolation** `${}` is THE most important syntax
2. **Closures** are like Python lambdas (you'll see `{ it * 2 }`)
3. **Maps** are like Python dicts (use dot notation)
4. **Lists** work like Python lists
5. **You don't need to be a Groovy expert** — just fluent

---

*This card summarizes Day 1. Review it before starting Day 2!*

**Progress: 1/28 days complete! 🎉**
