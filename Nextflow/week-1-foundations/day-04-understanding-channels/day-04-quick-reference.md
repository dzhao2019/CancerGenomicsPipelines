# Day 4: Channels and Workflows - Quick Reference

**Print this for quick lookup while building workflows!**

---

## 📋 Creating Channels (Foundation)

```groovy
// From files
Channel.fromPath("data/*.fastq")

// From values
Channel.from("s1", "s2", "s3")

// Single value (shared with all)
Channel.value("reference.fa")

// Paired files
Channel.fromFilePairs("reads/*_{1,2}.fastq.gz")

// From CSV
Channel.fromPath("metadata.csv").splitCsv(header:true)
```

---

## 🔄 Transforming Channels

```groovy
// Transform each item
channel.map { item -> transform(item) }

// Keep matching items
channel.filter { item -> condition(item) }

// Gather all items into one
channel.collect()

// Combine two channels
channel1.join(channel2)

// Split into branches
channel.branch {
    branch1: condition1
    branch2: condition2
}

// Chain operations
channel
    .filter { it.size() > 100 }
    .map { it.baseName }
    .collect()
```

---

## 🏗️ Workflow Template

```groovy
workflow {
    // 1. Create channels
    input_data = Channel.fromPath("data/*.fastq")
    reference = Channel.value("reference.fa")
    
    // 2. Transform (optional)
    data = input_data.map { f -> [f.baseName, f] }
    
    // 3. Pass to processes
    results1 = PROCESS1(data)
    results2 = PROCESS2(results1, reference)
    
    // 4. View or collect results
    results2.view()
}
```

---

## 💡 Common Patterns

### Pattern: Multiple Inputs

```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    reference = Channel.value("ref.fa")
    
    // reference used by all samples
    PROCESS(samples, reference)
}
```

### Pattern: Paired Data

```groovy
workflow {
    reads = Channel.fromFilePairs("*_{1,2}.fastq.gz")
    PROCESS(reads)
    // Each pair goes as unit
}
```

### Pattern: Sequential Processing

```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    
    step1 = PROCESS1(samples)
    step2 = PROCESS2(step1)
    step3 = PROCESS3(step2)
    
    step3.view()
}
```

### Pattern: Parallel Processes

```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    
    // Both run on same input
    qc = QC_PROCESS(samples)
    aligned = ALIGN_PROCESS(samples)
    
    aligned.view()
}
```

### Pattern: Collect All Items

```groovy
workflow {
    files = Channel.fromPath("*.bam")
    
    // Collect all BAM files
    all_bams = files.collect()
    
    // Pass all to merge process
    MERGE_PROCESS(all_bams)
}
```

---

## 📝 Channel Type Reference

| Type | Creation | Items | Usage |
|------|----------|-------|-------|
| Queue | `fromPath()` | Multiple | One at a time (parallelized) |
| Queue | `from()` | Multiple | One at a time (parallelized) |
| Value | `value()` | Single | Broadcast to all |
| Queue | `fromFilePairs()` | Multiple tuples | Paired data |

---

## 🔄 Transform Reference

```groovy
// Example channel: [1, 2, 3, 4, 5]

.map { it * 2 }
// Output: [2, 4, 6, 8, 10]

.filter { it > 3 }
// Output: [4, 5] (after map: [8, 10])

.collect()
// Output: [[2, 4, 6, 8, 10]] (single item)

.size()
// Output: 5 (number of items)
```

---

## 🎯 Connecting Processes

### Simple Connection
```groovy
workflow {
    input = Channel.fromPath("*.fastq")
    output = PROCESS(input)
}
```

### With Reference
```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    reference = Channel.value("ref.fa")
    
    output = PROCESS(samples, reference)
}
```

### Sequential Chain
```groovy
workflow {
    data = Channel.fromPath("*.fastq")
    
    out1 = PROC1(data)
    out2 = PROC2(out1)
    out3 = PROC3(out2)
}
```

### With Tuples
```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
        .map { f -> [f.baseName, f] }
    
    // Process receives: [sample_id, file]
    PROCESS(samples)
}
```

---

## 🐛 Common Issues & Fixes

| Problem | Cause | Fix |
|---------|-------|-----|
| No output | Process never called | Add `PROCESS(channel)` in workflow |
| Wrong data type | Using `path` instead of `val` | Check input types match |
| Files not found | Pattern doesn't match | Check path and pattern |
| Too many items | Using `.collect()` incorrectly | Only collect if you need all together |
| Process runs wrong times | Channel issue | Debug with `.view()` |

---

## 🔍 Debugging Techniques

```groovy
// See what's in a channel
channel.view()

// See with label
channel.view { "Item: $it" }

// Count items
channel.count().view()

// See file properties
files.view { f -> "File: ${f.name} Size: ${f.size()}" }

// Check if channel is empty
channel.ifEmpty { error("No items!") }
```

---

## 📊 Data Flow Timeline

**Workflow with 3 samples:**

```
Input files (3)
    ↓
Channel emits 3 items
    ↓
Process 1 (runs 3 times in parallel if resources allow)
    ↓
3 outputs
    ↓
Process 2 (runs 3 times in parallel)
    ↓
Final outputs
```

**With .collect():**

```
Process outputs (3 items)
    ↓
.collect() gathers all 3
    ↓
Process 3 (runs ONCE with all items)
    ↓
Single output
```

---

## ✅ Workflow Checklist

Before running your workflow:

- [ ] All input channels created
- [ ] All processes defined
- [ ] All processes called (not just defined)
- [ ] Output of one process used by next
- [ ] Reference channels as value if shared
- [ ] Sample IDs preserved in tuples
- [ ] .collect() only used when necessary

---

## 🎯 Most Important Rules

1. **Channels = Data Streams**
   - Multiple items = automatic parallelization
   - Single value = broadcast to all

2. **Pass Channels to Processes**
   - `channel.view()` doesn't run processes
   - `PROCESS(channel)` actually runs

3. **Preserve Metadata**
   - Use tuples to keep sample IDs
   - Sample ID should travel through pipeline

4. **Understand .collect()**
   - Without it: parallel (one item at a time)
   - With it: sequential (all items together)

5. **Test with .view()**
   - Debug channels before processes
   - See what's actually flowing

---

## 📋 Workflow Template (Copy & Modify)

```groovy
process STEP1 {
    input: tuple val(id), path(file)
    output: tuple val(id), path("*.output")
    script: """
        command ${file} > output.txt
    """
}

process STEP2 {
    input: tuple val(id), path(file)
    output: tuple val(id), path("*.final")
    script: """
        command2 ${file} > final.txt
    """
}

workflow {
    // Create channels
    input = Channel.fromPath("data/*.input")
        .map { f -> [f.baseName, f] }
    
    // Connect processes
    step1_out = STEP1(input)
    step2_out = STEP2(step1_out)
    
    // View results
    step2_out.view()
}
```

---

## 🚀 Quick Start: Three Workflows

### Simple (Single Process)
```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    PROCESS(samples)
}
```

### Chained (Two Processes)
```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    out1 = PROC1(samples)
    PROC2(out1)
}
```

### Complex (Multiple Paths)
```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    reference = Channel.value("ref.fa")
    
    qc = QC(samples)
    aligned = ALIGN(samples, reference)
    variants = CALL(aligned, reference)
    
    variants.view()
}
```

---

*Print this page and keep it handy while building workflows!*

**Next:** Day 5 - Subworkflows and best practices

