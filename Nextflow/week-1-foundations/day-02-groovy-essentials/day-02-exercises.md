# Day 2: Groovy Essentials - Hands-On Exercises

**Time**: 15 minutes for exercises + 10 minutes for solutions  
**Difficulty**: Easy to Moderate  
**Goal**: Practice Groovy concepts in realistic Nextflow contexts

---

## How to Use These Exercises

1. **Try first** - Attempt each exercise before looking at the solution
2. **Check your answer** - Compare with the solution provided
3. **Understand why** - Read the explanation to understand the concept
4. **Try variations** - Modify the exercise to test your understanding

---

## Exercise 1: String Interpolation (5 minutes)

### Part A: Basic Interpolation

**What will this print?**

```groovy
sample = "tumor_sample"
fastq = "sample_001.fastq"
result = "Processing ${sample} from file ${fastq}"
println result
```

**Your answer:**
```
[Write what you think prints here]
```

<details>
<summary>✓ Solution A</summary>

```
Processing tumor_sample from file sample_001.fastq
```

**Explanation:**
- `${sample}` is replaced with the value `"tumor_sample"`
- `${fastq}` is replaced with the value `"sample_001.fastq"`
- The resulting string is printed

**Key concept:** `${}` allows variable substitution in double-quoted strings.

</details>

---

### Part B: Interpolation with Expressions

**What will this print?**

```groovy
count = 1000000
processed = 750000
percent = (processed / count) * 100

message = "Processed ${processed} of ${count} reads (${percent}%)"
println message
```

**Your answer:**
```
[Write what you think prints here]
```

<details>
<summary>✓ Solution B</summary>

```
Processed 750000 of 1000000 reads (75.0%)
```

**Explanation:**
- `${processed}` → `750000`
- `${count}` → `1000000`
- `${percent}` → `75.0` (math is evaluated)
- Expressions inside `${}` are calculated

**Key concept:** You can put calculations and method calls inside `${}`.

</details>

---

### Part C: Interpolation with Methods

**What will this print?**

```groovy
sample_file = "sample_001.fastq"
message = "File name: ${sample_file.toUpperCase()}"
println message
```

**Your answer:**
```
[Write what you think prints here]
```

<details>
<summary>✓ Solution C</summary>

```
File name: SAMPLE_001.FASTQ
```

**Explanation:**
- `${sample_file.toUpperCase()}` calls the method on the string
- `.toUpperCase()` returns the uppercase version
- The result is interpolated

**Key concept:** Methods can be called inside `${}`.

</details>

---

## Exercise 2: Lists and Closures (5 minutes)

### Part A: Collecting (Transforming Lists)

**What will this print?**

```groovy
samples = ["sample1", "sample2", "sample3"]
fastq_files = samples.collect { sample -> "${sample}.fastq" }
println fastq_files
```

**Your answer:**
```
[Write what you think prints here]
```

<details>
<summary>✓ Solution A</summary>

```
[sample1.fastq, sample2.fastq, sample3.fastq]
```

**Explanation:**
- `.collect { sample -> ... }` transforms each item
- For each `sample`, it runs the closure
- The closure returns `"${sample}.fastq"`
- Results are collected into a new list

**Key concept:** `.collect()` with a closure is like Python's `map()` function.

**Python equivalent:**
```python
samples = ["sample1", "sample2", "sample3"]
fastq_files = [f"{sample}.fastq" for sample in samples]
# or: list(map(lambda sample: f"{sample}.fastq", samples))
```

</details>

---

### Part B: Filtering Lists

**What will this print?**

```groovy
numbers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
evens = numbers.findAll { number -> number % 2 == 0 }
println evens
```

**Your answer:**
```
[Write what you think prints here]
```

<details>
<summary>✓ Solution B</summary>

```
[2, 4, 6, 8, 10]
```

**Explanation:**
- `.findAll { ... }` keeps items where closure returns true
- `number % 2 == 0` is true only for even numbers
- Only matching items are kept in the result

**Key concept:** `.findAll()` with a closure is like Python's `filter()` or list comprehension with `if`.

**Python equivalent:**
```python
numbers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
evens = [n for n in numbers if n % 2 == 0]
# or: list(filter(lambda n: n % 2 == 0, numbers))
```

</details>

---

### Part C: Using `it` in Closures

**What will this print?**

```groovy
tools = ["fastqc", "bwa", "samtools"]
uppercase = tools.collect { it.toUpperCase() }
println uppercase
```

**Your answer:**
```
[Write what you think prints here]
```

<details>
<summary>✓ Solution C</summary>

```
[FASTQC, BWA, SAMTOOLS]
```

**Explanation:**
- When a closure has one parameter, you can use `it` instead of naming it
- `{ it.toUpperCase() }` is shorthand for `{ item -> item.toUpperCase() }`
- `it` refers to the current item

**Key concept:** `it` is the default parameter name in Groovy closures.

</details>

---

## Exercise 3: Maps (5 minutes)

### Part A: Creating and Accessing Maps

**What will this print?**

```groovy
sample = [id: "sample1", fastq: "s1.fastq", condition: "control"]

id = sample.id
fastq = sample["fastq"]
condition = sample.condition

println "ID: ${id}, File: ${fastq}, Condition: ${condition}"
```

**Your answer:**
```
[Write what you think prints here]
```

<details>
<summary>✓ Solution A</summary>

```
ID: sample1, File: s1.fastq, Condition: control
```

**Explanation:**
- `sample.id` accesses the value using dot notation
- `sample["fastq"]` accesses the value using bracket notation
- Both work the same way
- All values are interpolated into the string

**Key concept:** Maps are like Python dicts. Use either `map.key` or `map["key"]`.

</details>

---

### Part B: Iterating Over Maps

**What will this print?**

```groovy
config = [threads: 8, memory: "16GB", timeout: "2h"]

config.each { key, value ->
    println "${key}: ${value}"
}
```

**Your answer:**
```
[Write what you think prints here]
```

<details>
<summary>✓ Solution B</summary>

```
threads: 8
memory: 16GB
timeout: 2h
```

**Explanation:**
- `.each { key, value -> ... }` iterates over map entries
- For each entry, `key` is the key and `value` is the value
- The closure is executed once per entry

**Key concept:** Use `.each { key, value -> ... }` to iterate over maps.

**Python equivalent:**
```python
config = {"threads": 8, "memory": "16GB", "timeout": "2h"}
for key, value in config.items():
    print(f"{key}: {value}")
```

</details>

---

### Part C: Transforming Maps with Closures

**What will this print?**

```groovy
samples = [
    [id: "s1", reads: 1000000],
    [id: "s2", reads: 2500000],
    [id: "s3", reads: 1500000]
]

ids = samples.collect { sample -> sample.id }
println ids
```

**Your answer:**
```
[Write what you think prints here]
```

<details>
<summary>✓ Solution C</summary>

```
[s1, s2, s3]
```

**Explanation:**
- `samples` is a list of maps
- `.collect { sample -> sample.id }` extracts the `id` from each map
- Returns a list of just the IDs

**Key concept:** You can use `.collect()` on a list of maps to extract specific fields.

**Real Nextflow example:**
```groovy
// Extract sample IDs from metadata
sample_ids = metadata.collect { sample -> sample.id }
```

</details>

---

## Exercise 4: Reading Nextflow Code (5 minutes)

### Part A: Simple Process Reading

**What does this process do?**

```groovy
process ALIGN_READS {
    input:
        path fastq
        path reference
    
    output:
        path "*.bam"
    
    script:
        """
        bwa mem -t 4 ${reference} ${fastq} | samtools view -b > output.bam
        """
}
```

**Your answer:**
- **What is `fastq`?** [Your answer]
- **What does `${reference}` do?** [Your answer]
- **What is the output?** [Your answer]

<details>
<summary>✓ Solution A</summary>

**What is `fastq`?**
- It's an input file path. The Groovy variable `fastq` holds the path to the FASTQ file.

**What does `${reference}` do?**
- It interpolates the actual reference file path into the bash command.
- When this runs, `${reference}` is replaced with the actual file path (e.g., `/data/hg38.fa`).

**What is the output?**
- Any file ending in `.bam` (matched by the glob pattern `*.bam`)
- In this case, `output.bam`

**Key concept:** The Groovy variables are substituted into the bash script before it runs.

</details>

---

### Part B: Process with Parameters

**What will happen when this runs with `--quality 20`?**

```groovy
params {
    quality = 30
}

process FILTER_READS {
    input:
        path reads
    output:
        path "*.filtered.fastq"
    script:
        """
        seqtk seq -q ${params.quality} ${reads} > output.filtered.fastq
        """
}
```

**Your answer:**
- **What is `${params.quality}`?** [Your answer]
- **What command will run?** [Your answer]

<details>
<summary>✓ Solution B</summary>

**What is `${params.quality}`?**
- It's accessing the `quality` parameter from the `params` map
- With `--quality 20`, the value is 20
- With no parameter specified, the default is 30

**What command will run (with --quality 20)?**
```bash
seqtk seq -q 20 reads.fastq > output.filtered.fastq
```

**Explanation:** 
- `${params.quality}` is interpolated into the command
- The final command uses the quality threshold from parameters

**Key concept:** `params` is a special map in Nextflow. Access it with `${params.key}`.

</details>

---

### Part C: Channel Operations

**What does this do?**

```groovy
workflow {
    samples = Channel.fromPath("data/*.fastq")
    
    names = samples.map { file -> file.baseName }
    
    names.view()
}
```

**Your answer:**
- **What does `Channel.fromPath()`?** [Your answer]
- **What does `.map() { file -> file.baseName }`?** [Your answer]
- **What does `.view()`?** [Your answer]

<details>
<summary>✓ Solution C</summary>

**What does `Channel.fromPath()`?**
- Creates a channel with all files matching the pattern `data/*.fastq`
- If there are 3 files: `sample1.fastq`, `sample2.fastq`, `sample3.fastq`
- The channel will emit 3 items

**What does `.map() { file -> file.baseName }`?**
- Transforms each file into just its base name (without path and extension)
- `sample1.fastq` → `sample1`
- `sample2.fastq` → `sample2`
- `sample3.fastq` → `sample3`

**What does `.view()`?**
- Prints each item to the console for debugging
- Output would be:
  ```
  sample1
  sample2
  sample3
  ```

**Key concept:** Channels are streams of data. `.map()` transforms each item, `.view()` prints for debugging.

</details>

---

## Exercise 5: Debugging Groovy Errors (5 minutes)

### Part A: Spot the Error

**What's wrong with this code?**

```groovy
samples = ['sample1', 'sample2', 'sample3']
result = samples.collect { it.uppercase }
println result
```

**Your answer:**
```
[What's the error?]
```

<details>
<summary>✓ Solution A</summary>

**Error:** `it.uppercase` is missing parentheses!

**Wrong:**
```groovy
result = samples.collect { it.uppercase }  // Method needs ()
```

**Right:**
```groovy
result = samples.collect { it.toUpperCase() }  // Correct
```

Or even shorter:
```groovy
result = samples.collect { it.toUpperCase() }
// Result: [SAMPLE1, SAMPLE2, SAMPLE3]
```

**Explanation:** 
- In Groovy, methods must be called with parentheses
- `it.uppercase` tries to access a property, not call a method
- `.toUpperCase()` is the actual method name

**Key concept:** Methods need `()` even if they have no parameters.

</details>

---

### Part B: Spot the Error

**What's wrong with this code?**

```groovy
name = "Alice"
message = 'Hello ${name}'
println message
```

**Your answer:**
```
[What's the error or unexpected behavior?]
```

<details>
<summary>✓ Solution B</summary>

**Issue:** Single quotes don't interpolate!

**Output:**
```
Hello ${name}
```

Not:
```
Hello Alice
```

**Explanation:** 
- Single quotes `'...'` are literal strings in Groovy (no interpolation)
- Double quotes `"..."` allow interpolation with `${}`
- The `${}` is printed literally, not evaluated

**Fix:**
```groovy
message = "Hello ${name}"  // Use double quotes
println message             // Prints: Hello Alice
```

**Key concept:** Use double quotes `"..."` for interpolation, single quotes `'...'` for literals.

</details>

---

### Part C: Spot the Error

**What's wrong with this code?**

```groovy
config = {"name": "Alice", "age": 30}
println config.name
```

**Your answer:**
```
[What's the error?]
```

<details>
<summary>✓ Solution C</summary>

**Error:** Wrong map syntax! This is Python dict syntax, not Groovy.

**Wrong:**
```groovy
config = {"name": "Alice", "age": 30}  // Python dict syntax
```

**Right:**
```groovy
config = [name: "Alice", age: 30]      // Groovy map syntax
println config.name                     // Prints: Alice
```

**Explanation:**
- Groovy maps use `[key: value]` syntax
- Python dicts use `{"key": value}` syntax
- They look similar but are different

**Key concept:** Groovy maps: `[key: value]`
Python dicts: `{"key": value}`

</details>

---

## Exercise 6: Writing Simple Groovy (5 minutes)

### Part A: Write a Closure

**Write a closure that filters numbers > 5**

```groovy
numbers = [1, 3, 5, 7, 9, 11]

# Write your closure here
filtered = numbers.findAll { ??? }

println filtered
```

**Your answer:**
```groovy
filtered = numbers.findAll { ??? }
```

<details>
<summary>✓ Solution A</summary>

```groovy
filtered = numbers.findAll { it > 5 }
```

**Output:**
```
[7, 9, 11]
```

**Explanation:**
- `{ it > 5 }` is the closure
- It returns true for numbers > 5
- `.findAll()` keeps only matching items

**Alternative (more explicit):**
```groovy
filtered = numbers.findAll { number -> number > 5 }
```

Both work the same way.

</details>

---

### Part B: Interpolate Values

**Complete the string interpolation**

```groovy
sample_id = "tumor_001"
reads = 5000000
message = "???"  # Write the string

println message
# Should print: Sample tumor_001 has 5000000 reads
```

**Your answer:**
```groovy
message = ???
```

<details>
<summary>✓ Solution B</summary>

```groovy
message = "Sample ${sample_id} has ${reads} reads"
```

**Output:**
```
Sample tumor_001 has 5000000 reads
```

**Explanation:**
- Use double quotes for interpolation
- `${sample_id}` becomes `"tumor_001"`
- `${reads}` becomes `5000000`

</details>

---

### Part C: Create and Use a Map

**Create a map for process metadata and access values**

```groovy
# Create a map with: id="s1", fastq="s1.fastq", cores=4

# Access the id and cores

# Print: "Process s1 on 4 cores"
```

**Your answer:**
```groovy
[Write your code here]
```

<details>
<summary>✓ Solution C</summary>

```groovy
process_info = [id: "s1", fastq: "s1.fastq", cores: 4]

id = process_info.id
cores = process_info.cores

println "Process ${id} on ${cores} cores"
```

**Output:**
```
Process s1 on 4 cores
```

**Or more concisely:**
```groovy
process_info = [id: "s1", fastq: "s1.fastq", cores: 4]
println "Process ${process_info.id} on ${process_info.cores} cores"
```

**Explanation:**
- Create map with `[key: value]` syntax
- Access with `map.key` or `map["key"]`
- Interpolate with `${...}`

</details>

---

## Exercise 7: Real Nextflow Pattern (10 minutes)

### Challenge: Understand a Real Process

**Given this process, answer the questions:**

```groovy
params {
    output_dir = "results"
}

process QUALITY_REPORT {
    publishDir "${params.output_dir}/qc"
    
    input:
        tuple val(sample_id), path(fastq), val(read_count)
    
    output:
        path "*.report"
    
    script:
        """
        echo "Sample: ${sample_id}" > report.txt
        echo "File: ${fastq}" >> report.txt
        echo "Reads: ${read_count}" >> report.txt
        echo "Quality report for ${sample_id}" > ${sample_id}.report
        cat report.txt >> ${sample_id}.report
        """
}
```

**Questions:**

1. **What type of data is in the `input: tuple`?**
   - Answer: [Your answer]

2. **What does `publishDir "${params.output_dir}/qc"` do?**
   - Answer: [Your answer]

3. **What files will be created?**
   - Answer: [Your answer]

4. **What is `${sample_id}` replaced with?**
   - Answer: [Your answer]

<details>
<summary>✓ Solutions</summary>

**1. What type of data is in the `input: tuple`?**
- `val(sample_id)` - A string value with the sample ID
- `path(fastq)` - A file path to the FASTQ file
- `val(read_count)` - A number value with the read count

**Explanation:** `tuple` groups related data together. Each item is labeled with its type (val or path).

**2. What does `publishDir "${params.output_dir}/qc"` do?**
- Copies output files to `results/qc/` directory
- Uses the `output_dir` parameter
- With default params, outputs go to `results/qc/`

**Explanation:** `publishDir` is a Nextflow directive that specifies where to copy outputs.

**3. What files will be created?**
- `${sample_id}.report` - The main report file
- `report.txt` - Created during execution (but not published unless in output)

**Output:** The `*.report` pattern means only files matching `.report` are published.

**4. What is `${sample_id}` replaced with?**
- The actual sample ID passed to the process
- Example: If `sample_id = "tumor_1"`, then `${sample_id}` becomes `tumor_1`
- The output file would be `tumor_1.report`

**Explanation:** String interpolation replaces `${}` with the actual value.

---

**Real-world example:**
```groovy
// When called with:
// sample_id = "sample_001"
// fastq = "/path/to/sample_001.fastq"
// read_count = 5000000

// The script becomes:
"""
echo "Sample: sample_001" > report.txt
echo "File: /path/to/sample_001.fastq" >> report.txt
echo "Reads: 5000000" >> report.txt
echo "Quality report for sample_001" > sample_001.report
cat report.txt >> sample_001.report
"""
```

</details>

---

## Summary: What You've Practiced

✅ String interpolation (the most important!)  
✅ Lists and closures (transform and filter)  
✅ Maps (store and access data)  
✅ Reading Nextflow code  
✅ Spotting and fixing errors  
✅ Writing simple Groovy  
✅ Understanding real processes  

---

## 🎓 Self-Check

**Can you do these without looking at solutions?**

- [ ] Write a closure that doubles numbers
- [ ] Interpolate three variables into a string
- [ ] Access values from a map using both `.` and `["key"]` notation
- [ ] Explain what `.collect()` does
- [ ] Spot the difference between `'...'` and `"..."`
- [ ] Understand a simple Nextflow process
- [ ] Fix a Groovy syntax error

If you can do most of these, you're ready for Day 3!

---

## 🚀 Ready for Day 3?

You now understand Groovy well enough to:
- Read any Nextflow workflow
- Understand what's happening
- Write simple processes
- Debug common errors

Tomorrow you'll write your first Nextflow process. Everything you learned today will make it straightforward!

---

*Day 2 of 28 complete! You're building the skills to write production-ready workflows.*
