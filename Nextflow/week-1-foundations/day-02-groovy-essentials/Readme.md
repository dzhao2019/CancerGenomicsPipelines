# Day 2: Groovy Essentials for Nextflow

**Learning Time**: 30 minutes  
**Prerequisites**: Day 1 completed, Basic Python knowledge  
**Goal**: Read and write Groovy code comfortably for Nextflow workflows

---

## 📖 Introduction (3 minutes)

Welcome to Day 2! Today you'll learn Groovy—but don't worry, you're not becoming a Groovy expert. You're learning just enough to read and write Nextflow workflows comfortably.

**The Good News**: If you know Python, you already understand 80% of what you need. Groovy is like Python's Java-based cousin—similar ideas, slightly different syntax.

### What You'll Learn Today

- Core Groovy syntax you'll use in Nextflow
- String interpolation (the `${}` pattern you'll see everywhere)
- Closures (Groovy's version of lambda functions)
- Lists and Maps (similar to Python lists and dicts)
- Key differences from Python that matter
- What you DON'T need to learn (most of Groovy!)

### What You WON'T Learn

You're not learning:
- ❌ Object-oriented Groovy programming
- ❌ Advanced Groovy features
- ❌ Groovy frameworks
- ❌ Java interoperability details

You're learning:
- ✅ Reading Nextflow workflow syntax
- ✅ Writing simple Groovy expressions
- ✅ Working with data structures in Nextflow
- ✅ Understanding patterns you'll see daily

---

## 🎯 The Big Picture: Why Groovy? (2 minutes)

**Why not just use Python for Nextflow?**

Nextflow chose Groovy because:
1. **JVM-based**: Runs on the Java Virtual Machine (portable, robust)
2. **Dynamic**: Similar feel to Python (no type declarations needed)
3. **Concise**: Less verbose than Java
4. **DSL-friendly**: Great for creating Domain Specific Languages (which Nextflow is)

**For you**: It means learning some new syntax, but the logic is familiar.

---

## 🔑 Key Concepts with Examples (15 minutes)

### 1. Variables and Basic Types (2 minutes)

**Python**:
```python
name = "Alice"
age = 30
height = 5.6
is_scientist = True
```

**Groovy**:
```groovy
name = "Alice"
age = 30
height = 5.6
is_scientist = true  // lowercase true/false
```

**Key Differences**:
- Groovy uses `true/false` (lowercase) instead of `True/False`
- Groovy uses `null` instead of `None`
- Semicolons are optional (usually omitted)
- Type declarations are optional

```groovy
// All of these work in Groovy
def name = "Alice"        // 'def' means dynamic type
String name = "Alice"     // Explicit type (optional)
name = "Alice"            // Type inferred (most common in Nextflow)
```

**In Nextflow, you'll see**:
```groovy
params.input = "data/"
params.output = "results/"
params.threads = 8
```

### 2. String Interpolation - THE MOST IMPORTANT PATTERN (3 minutes)

This is crucial because you'll use it constantly in Nextflow!

**Python**:
```python
name = "Alice"
age = 30

# F-strings (modern Python)
message = f"Hello, {name}! You are {age} years old."

# Format method
message = "Hello, {}! You are {} years old.".format(name, age)

# Old style
message = "Hello, %s! You are %d years old." % (name, age)
```

**Groovy**:
```groovy
name = "Alice"
age = 30

// String interpolation with ${}
message = "Hello, ${name}! You are ${age} years old."

// For simple variables, you can skip the braces
message = "Hello, $name! You are $age years old."
```

**When to use braces**:
```groovy
// ALWAYS use braces for:
// 1. Expressions
output = "Result: ${x + y}"

// 2. Property access
output = "File: ${file.name}"

// 3. When variable name is followed by more text
output = "sample_${id}_processed.bam"  // Correct
output = "sample_$id_processed.bam"    // Wrong! Looks for variable "id_processed"
```

**In Nextflow scripts** (you'll see this everywhere):
```groovy
process align {
    input:
    val sample_id
    path reads
    
    output:
    path "${sample_id}.bam"
    
    script:
    """
    bwa mem reference.fa ${reads} > ${sample_id}.bam
    samtools sort ${sample_id}.bam -o ${sample_id}.sorted.bam
    """
}
```

**Pro Tip**: In Nextflow `script` blocks (triple quotes), you use `${variable}` to reference Groovy variables. This is THE pattern you'll use most!

### 3. Lists (Like Python Lists) (2 minutes)

**Python**:
```python
# Creating lists
samples = ["sample1", "sample2", "sample3"]
numbers = [1, 2, 3, 4, 5]

# Accessing elements
first = samples[0]
last = samples[-1]

# List operations
samples.append("sample4")
length = len(samples)

# List comprehension
doubled = [x * 2 for x in numbers]
```

**Groovy**:
```groovy
// Creating lists
samples = ["sample1", "sample2", "sample3"]
numbers = [1, 2, 3, 4, 5]

// Accessing elements
first = samples[0]
last = samples[-1]  // Yes, negative indexing works!

// List operations
samples << "sample4"  // << is like append
samples.add("sample5")  // .add() also works
length = samples.size()

// Using collect (like map)
doubled = numbers.collect { it * 2 }
```

**Key Differences**:
- `<<` operator adds to list (common in Groovy)
- `size()` instead of `len()`
- `collect {}` is like Python's list comprehension

**Common Groovy List Methods**:
```groovy
samples = ["sample1", "sample2", "sample3"]

// Filter
long_names = samples.findAll { it.length() > 7 }

// Map/transform
uppercase = samples.collect { it.toUpperCase() }

// Each (like for loop)
samples.each { println it }

// Join
text = samples.join(", ")  // "sample1, sample2, sample3"
```

**In Nextflow**:
```groovy
// You'll often see
Channel.from(["sample1", "sample2", "sample3"])

// Or using ranges
Channel.from(1..10)  // Numbers 1 through 10
```

### 4. Maps (Like Python Dictionaries) (2 minutes)

**Python**:
```python
# Creating dictionaries
config = {
    "input": "data/",
    "output": "results/",
    "threads": 8
}

# Accessing values
input_dir = config["input"]
threads = config.get("threads", 4)

# Adding/updating
config["genome"] = "hg38"
```

**Groovy**:
```groovy
// Creating maps
config = [
    input: "data/",
    output: "results/",
    threads: 8
]

// Accessing values (multiple ways)
input_dir = config["input"]     // Bracket notation
input_dir = config.input        // Dot notation
threads = config.get("threads", 4)

// Adding/updating
config.genome = "hg38"
config["reference"] = "ref.fa"
```

**Key Differences**:
- No quotes needed for keys (unless they have spaces/special chars)
- Can use dot notation for access
- Square brackets `[]` for maps, not curly braces `{}`

**In Nextflow Parameters**:
```groovy
// Very common pattern
params.input = "data/*.fastq"
params.output = "results/"
params.genome = "hg38"
params.threads = 8

// Accessing in processes
process align {
    cpus params.threads
    
    script:
    """
    bwa mem -t ${params.threads} reference.fa reads.fastq
    """
}
```

### 5. Closures - Groovy's Lambda Functions (3 minutes)

This is where Groovy differs most from Python, but it's actually quite elegant!

**Python Lambda**:
```python
# Simple lambda
double = lambda x: x * 2

# With map
numbers = [1, 2, 3, 4, 5]
doubled = list(map(lambda x: x * 2, numbers))

# With filter
evens = list(filter(lambda x: x % 2 == 0, numbers))
```

**Groovy Closures**:
```groovy
// Closure (like lambda)
double = { x -> x * 2 }

// With collect (like map)
numbers = [1, 2, 3, 4, 5]
doubled = numbers.collect { x -> x * 2 }

// Implicit 'it' parameter (super common!)
doubled = numbers.collect { it * 2 }

// With findAll (like filter)
evens = numbers.findAll { it % 2 == 0 }
```

**The `it` keyword**: When a closure has one parameter, you can use `it` instead of naming it.

```groovy
// These are equivalent:
samples.collect { sample -> sample.toUpperCase() }
samples.collect { it.toUpperCase() }

// These are equivalent:
numbers.findAll { n -> n > 5 }
numbers.findAll { it > 5 }
```

**Multi-line Closures**:
```groovy
// For complex operations
results = samples.collect { sample ->
    def trimmed = sample.trim()
    def upper = trimmed.toUpperCase()
    return "${upper}_processed"
}

// Last expression is automatically returned
results = samples.collect { sample ->
    def trimmed = sample.trim()
    trimmed.toUpperCase()  // Returned automatically
}
```

**In Nextflow** (you'll see closures everywhere):
```groovy
// Transforming channel data
Channel
    .fromPath("*.fastq")
    .map { file -> [file.baseName, file] }
    .view { name, file -> "Processing: ${name}" }

// Using 'it' for simple cases
Channel
    .fromPath("*.fastq")
    .map { it.baseName }
    .view { "Sample: ${it}" }
```

### 6. Strings: Single vs Double Quotes (2 minutes)

**Python**: Single and double quotes are mostly interchangeable
```python
name = 'Alice'  # Same as "Alice"
message = f"Hello, {name}"  # F-string for interpolation
```

**Groovy**: Quote type matters!
```groovy
// Single quotes: NO interpolation (literal string)
message = 'Hello, ${name}'  // Prints: Hello, ${name}

// Double quotes: YES interpolation
message = "Hello, ${name}"  // Prints: Hello, Alice

// Triple quotes: Multi-line with interpolation
script = """
    echo "Processing ${sample}"
    bwa mem ${reference} ${reads}
    """
```

**In Nextflow**, you'll see:
```groovy
process example {
    script:
    """
    # This is a multi-line shell script
    # ${variables} are interpolated from Groovy
    echo "Sample: ${sample_id}"
    fastqc ${input_file}
    """
}
```

**Important Rule**: 
- `"""triple quotes"""` for multi-line scripts (Nextflow script blocks)
- `"double quotes"` for single-line strings with interpolation
- `'single quotes'` for literal strings (rare in Nextflow)

### 7. Control Flow - Familiar Territory (2 minutes)

**If/Else**:

Python:
```python
if score > 90:
    grade = "A"
elif score > 80:
    grade = "B"
else:
    grade = "C"
```

Groovy:
```groovy
if (score > 90) {
    grade = "A"
} else if (score > 80) {
    grade = "B"
} else {
    grade = "C"
}
```

**For Loops**:

Python:
```python
for sample in samples:
    print(sample)

for i in range(10):
    print(i)
```

Groovy:
```groovy
// For-each style
for (sample in samples) {
    println sample
}

// Using each (more common in Groovy)
samples.each { sample ->
    println sample
}

// Ranges
for (i in 0..9) {
    println i
}

(0..9).each { println it }
```

**In Nextflow**, you typically use channels and operators instead of explicit loops, but you'll see these patterns in script blocks:

```groovy
process example {
    script:
    """
    # Shell script can have normal bash loops
    for file in *.bam; do
        echo \$file  # Note: \$ for bash variables
    done
    
    # But Groovy ${} for Nextflow variables
    echo "Sample ID: ${sample_id}"
    """
}
```

### 8. What You DON'T Need to Know (1 minute)

**Relax—you can skip**:
- Classes and object-oriented programming
- Groovy's advanced metaprogramming
- Gradle and Groovy build tools
- Most of the Groovy standard library
- Type systems and generics
- Exception handling details

**Focus on**:
- String interpolation `${}`
- Lists and Maps
- Closures with `it`
- Reading Nextflow examples

---

## 🔗 Python Connection: Side-by-Side Comparison (3 minutes)

Let's see a complete example in both languages:

### Task: Process a list of samples

**Python Version**:
```python
samples = ["sample1", "sample2", "sample3"]

# Transform: add suffix
processed = [f"{s}_processed.bam" for s in samples]

# Filter: only long names
long_names = [s for s in samples if len(s) > 7]

# Print each
for sample in samples:
    print(f"Processing {sample}")

# Dictionary of sample info
sample_info = {
    "sample1": {"type": "tumor", "batch": 1},
    "sample2": {"type": "normal", "batch": 1},
    "sample3": {"type": "tumor", "batch": 2}
}

# Access nested data
for sample, info in sample_info.items():
    print(f"{sample}: {info['type']} from batch {info['batch']}")
```

**Groovy Version**:
```groovy
samples = ["sample1", "sample2", "sample3"]

// Transform: add suffix
processed = samples.collect { "${it}_processed.bam" }

// Filter: only long names
long_names = samples.findAll { it.length() > 7 }

// Print each
samples.each { sample ->
    println "Processing ${sample}"
}

// Map of sample info
sample_info = [
    sample1: [type: "tumor", batch: 1],
    sample2: [type: "normal", batch: 1],
    sample3: [type: "tumor", batch: 2]
]

// Access nested data
sample_info.each { sample, info ->
    println "${sample}: ${info.type} from batch ${info.batch}"
}
```

### Common Patterns Translation Table

| Task | Python | Groovy |
|------|--------|--------|
| String interpolation | `f"{var}"` | `"${var}"` or `"$var"` |
| List comprehension | `[x*2 for x in list]` | `list.collect { it*2 }` |
| Filter list | `[x for x in list if x>5]` | `list.findAll { it>5 }` |
| For each | `for x in list: print(x)` | `list.each { println it }` |
| Dictionary | `{"key": "value"}` | `[key: "value"]` |
| Access dict | `dict["key"]` or `dict.get("key")` | `dict["key"]` or `dict.key` |
| Length | `len(list)` | `list.size()` |
| Append to list | `list.append(x)` | `list << x` or `list.add(x)` |
| Join strings | `", ".join(list)` | `list.join(", ")` |
| True/False | `True/False` | `true/false` |
| None/null | `None` | `null` |

---

## 💻 Hands-On Exercises (10 minutes)

### Exercise 1: String Interpolation Practice (3 minutes)

Convert these Python strings to Groovy with proper interpolation:

**Python**:
```python
sample_id = "patient001"
replicate = 2
file_type = "fastq"

# Task 1: Create filename
filename = f"{sample_id}_rep{replicate}.{file_type}"

# Task 2: Create command
threads = 8
command = f"fastqc -t {threads} {filename}"

# Task 3: Create output path
output_dir = "/results"
output_path = f"{output_dir}/{sample_id}/qc_rep{replicate}.html"
```

**Your Groovy Version**:
```groovy
sample_id = "patient001"
replicate = 2
file_type = "fastq"

// Task 1: Create filename

// Task 2: Create command

// Task 3: Create output path
```

<details>
<summary>Click to see solution</summary>

```groovy
sample_id = "patient001"
replicate = 2
file_type = "fastq"

// Task 1: Create filename
filename = "${sample_id}_rep${replicate}.${file_type}"
// Result: "patient001_rep2.fastq"

// Task 2: Create command
threads = 8
command = "fastqc -t ${threads} ${filename}"
// Result: "fastqc -t 8 patient001_rep2.fastq"

// Task 3: Create output path
output_dir = "/results"
output_path = "${output_dir}/${sample_id}/qc_rep${replicate}.html"
// Result: "/results/patient001/qc_rep2.html"

// Alternative for simple cases (no braces needed):
filename = "${sample_id}_rep${replicate}.$file_type"  // Works!
// But safer to always use braces for consistency
```

**Key Points**:
- Use `${}` for all variable interpolation
- Braces are required when variable name runs into other text
- Double quotes enable interpolation, single quotes don't
</details>

### Exercise 2: Lists and Closures (4 minutes)

Convert these Python operations to Groovy:

**Python**:
```python
samples = ["sample1", "sample2", "sample3", "sample4", "sample5"]

# Task 1: Create list of filenames like "sample1.bam"
bam_files = [f"{s}.bam" for s in samples]

# Task 2: Filter samples that contain "1" or "2"
early_samples = [s for s in samples if "1" in s or "2" in s]

# Task 3: Print each sample with index
for i, sample in enumerate(samples):
    print(f"Sample {i}: {sample}")

# Task 4: Create uppercase versions
upper_samples = [s.upper() for s in samples]

# Task 5: Get samples longer than 7 characters
long_names = [s for s in samples if len(s) > 7]
```

**Your Groovy Version**:
```groovy
samples = ["sample1", "sample2", "sample3", "sample4", "sample5"]

// Task 1: Create list of filenames like "sample1.bam"

// Task 2: Filter samples that contain "1" or "2"

// Task 3: Print each sample with index

// Task 4: Create uppercase versions

// Task 5: Get samples longer than 7 characters
```

<details>
<summary>Click to see solution</summary>

```groovy
samples = ["sample1", "sample2", "sample3", "sample4", "sample5"]

// Task 1: Create list of filenames like "sample1.bam"
bam_files = samples.collect { "${it}.bam" }
// Result: ["sample1.bam", "sample2.bam", "sample3.bam", ...]

// Alternative with explicit parameter name:
bam_files = samples.collect { sample -> "${sample}.bam" }

// Task 2: Filter samples that contain "1" or "2"
early_samples = samples.findAll { it.contains("1") || it.contains("2") }
// Result: ["sample1", "sample2"]

// Task 3: Print each sample with index
samples.eachWithIndex { sample, i ->
    println "Sample ${i}: ${sample}"
}
// Note: eachWithIndex gives you both element and index

// Task 4: Create uppercase versions
upper_samples = samples.collect { it.toUpperCase() }
// Result: ["SAMPLE1", "SAMPLE2", "SAMPLE3", ...]

// Task 5: Get samples longer than 7 characters
long_names = samples.findAll { it.length() > 7 }
// Result: [] (none are longer than 7 in this case)

// If we had longer names:
test = ["sample1", "sample2", "verylongsample"]
long_names = test.findAll { it.length() > 7 }
// Result: ["verylongsample"]
```

**Groovy Closure Patterns**:
- `collect { }` = transform each element (like Python map)
- `findAll { }` = filter elements (like Python filter)
- `each { }` = iterate (like Python for loop)
- `eachWithIndex { element, index -> }` = iterate with index
- Use `it` for simple single-parameter closures
- Use named parameters `{ x -> }` for clarity with multiple operations
</details>

### Exercise 3: Working with Maps (3 minutes)

Convert this Python dictionary manipulation to Groovy:

**Python**:
```python
# Sample metadata
samples = {
    "sample1": {
        "type": "tumor",
        "patient": "P001",
        "batch": 1
    },
    "sample2": {
        "type": "normal",
        "patient": "P001",
        "batch": 1
    },
    "sample3": {
        "type": "tumor",
        "patient": "P002",
        "batch": 2
    }
}

# Task 1: Access sample1's type
sample1_type = samples["sample1"]["type"]

# Task 2: Add a new sample
samples["sample4"] = {
    "type": "normal",
    "patient": "P002",
    "batch": 2
}

# Task 3: Print all tumor samples
for sample, info in samples.items():
    if info["type"] == "tumor":
        print(f"{sample} is from patient {info['patient']}")

# Task 4: Create list of all patients
patients = [info["patient"] for sample, info in samples.items()]
unique_patients = list(set(patients))
```

**Your Groovy Version**:
```groovy
// Sample metadata
samples = [
    sample1: [
        type: "tumor",
        patient: "P001",
        batch: 1
    ],
    sample2: [
        type: "normal",
        patient: "P001",
        batch: 1
    ],
    sample3: [
        type: "tumor",
        patient: "P002",
        batch: 2
    ]
]

// Task 1: Access sample1's type

// Task 2: Add a new sample

// Task 3: Print all tumor samples

// Task 4: Create list of all patients and get unique ones
```

<details>
<summary>Click to see solution</summary>

```groovy
// Sample metadata
samples = [
    sample1: [
        type: "tumor",
        patient: "P001",
        batch: 1
    ],
    sample2: [
        type: "normal",
        patient: "P001",
        batch: 1
    ],
    sample3: [
        type: "tumor",
        patient: "P002",
        batch: 2
    ]
]

// Task 1: Access sample1's type (multiple ways!)
sample1_type = samples.sample1.type           // Dot notation (cleanest)
sample1_type = samples["sample1"]["type"]     // Bracket notation
sample1_type = samples.sample1["type"]        // Mixed
// Result: "tumor"

// Task 2: Add a new sample
samples.sample4 = [
    type: "normal",
    patient: "P002",
    batch: 2
]
// Or using brackets:
samples["sample4"] = [type: "normal", patient: "P002", batch: 2]

// Task 3: Print all tumor samples
samples.each { sample, info ->
    if (info.type == "tumor") {
        println "${sample} is from patient ${info.patient}"
    }
}
// Output:
// sample1 is from patient P001
// sample3 is from patient P002

// Task 4: Create list of all patients and get unique ones
patients = samples.collect { sample, info -> info.patient }
// Result: ["P001", "P001", "P002", "P002"]

unique_patients = patients.unique()
// Result: ["P001", "P002"]

// Or in one line:
unique_patients = samples.collect { sample, info -> info.patient }.unique()
```

**Key Map Patterns**:
- Access with dots: `map.key.subkey`
- Access with brackets: `map["key"]["subkey"]`
- Iterate: `map.each { key, value -> }`
- Transform values: `map.collect { key, value -> }`
- No quotes needed for keys (unless special characters)
</details>

---

## 🤔 Reflection Activity (5 minutes)

### Question 1: Groovy vs Python Comfort Level

Look back at the exercises. Which Groovy concepts feel:

**Comfortable** (similar to Python):
- [ ] Variables and basic types
- [ ] Lists
- [ ] If/else statements
- [ ] For loops

**Slightly Different** (need practice):
- [ ] String interpolation with `${}`
- [ ] Maps (dictionaries)
- [ ] Dot notation for map access
- [ ] Method names (size vs len)

**New and Different** (need more practice):
- [ ] Closures and `it`
- [ ] `collect {}` and `findAll {}`
- [ ] The `<<` operator
- [ ] Triple-quoted strings

### Question 2: Pattern Practice

Write down these common patterns in your own words:

**String with variable**:
```groovy
// How would you write: "sample_{id}_processed.bam"?
// Your answer:
```

**Transform a list**:
```groovy
// How would you convert ["a", "b", "c"] to ["A", "B", "C"]?
// Your answer:
```

**Filter a list**:
```groovy
// How would you get numbers > 5 from [1,2,3,6,7,8]?
// Your answer:
```

<details>
<summary>Check your answers</summary>

```groovy
// String with variable
filename = "${sample_id}_processed.bam"

// Transform a list
uppercase = ["a", "b", "c"].collect { it.toUpperCase() }

// Filter a list
big_numbers = [1,2,3,6,7,8].findAll { it > 5 }
// Result: [6, 7, 8]
```
</details>

### Question 3: Real Nextflow Example

Here's a real Nextflow snippet. Can you understand it now?

```groovy
process trimReads {
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq")
    
    script:
    """
    trimmomatic SE ${reads} ${sample_id}_trimmed.fastq LEADING:3 TRAILING:3
    """
}

workflow {
    samples = Channel.fromFilePairs("data/*_{1,2}.fastq")
    trimmed = trimReads(samples)
}
```

**Questions to test your understanding**:
1. What does `${sample_id}` do in the output path?
2. What does the triple-quoted `"""..."""` section contain?
3. How is `${reads}` different from `$reads`?

<details>
<summary>Answers</summary>

1. **String interpolation**: Inserts the value of the `sample_id` variable into the filename. If `sample_id = "patient001"`, the output file is `"patient001_trimmed.fastq"`

2. **Shell script**: The triple-quoted section is a multi-line shell script that will be executed. It contains the actual trimmomatic command.

3. **Safety**: `${reads}` clearly marks where the variable starts and ends. While `$reads` might work in simple cases, `${reads}` is safer and more explicit, especially when the variable is followed by more text or when accessing properties like `${reads.name}`

**Key insight**: You can now read Nextflow processes! You understand:
- The `input:` declares what data comes in
- The `output:` declares what's produced
- The `script:` contains shell commands with Groovy interpolation
- `${variable}` injects Groovy values into the script
</details>

---

## 📝 Key Takeaways

Before moving to Day 3, make sure you understand:

✅ **String interpolation with `${variable}`** - You'll use this constantly  
✅ **Lists with `collect {}` and `findAll {}`** - Core data manipulation  
✅ **The `it` keyword** - Implicit parameter in closures  
✅ **Maps with dot and bracket notation** - For configuration and metadata  
✅ **Triple quotes `"""`** - For multi-line scripts in Nextflow  
✅ **Groovy is similar to Python** - Same logic, slightly different syntax  

### The Mental Model

Think of Groovy as Python's cousin who:
- Uses `${}` instead of `f"{}"`
- Says `collect` instead of list comprehension
- Says `findAll` instead of filter
- Likes both `map.key` and `map["key"]`
- Uses `it` as a shortcut in closures
- Is less picky about semicolons and parentheses

### What You Actually Need to Remember

**The Big 3 Patterns**:

1. **String interpolation**:
   ```groovy
   "${variable}"  // Always safe
   ```

2. **List transformation**:
   ```groovy
   list.collect { it.transform() }  // Like map
   list.findAll { it.condition }    // Like filter
   ```

3. **Map access**:
   ```groovy
   config.key        // Dot notation
   config["key"]     // Bracket notation
   ```

Everything else you can look up as needed!

---

## 🎯 Ready for Day 3?

Tomorrow, you'll write your **first Nextflow process**! You now have all the Groovy knowledge you need to understand process syntax.

### Quick Preview - You'll Understand This Tomorrow!

```groovy
process sayHello {
    input:
    val name
    
    output:
    path "greeting.txt"
    
    script:
    """
    echo "Hello, ${name}!" > greeting.txt
    """
}
```

You can already identify:
- `val name` - input parameter
- `${name}` - string interpolation
- `"""..."""` - multi-line script block
- `path "greeting.txt"` - output file

### Groovy Quick Reference Card

Save this for reference:

```groovy
// Variables
name = "Alice"
count = 42

// String interpolation
message = "Hello, ${name}"
filename = "${id}_processed.bam"

// Lists
samples = ["s1", "s2", "s3"]
samples << "s4"                    // Add item
transformed = samples.collect { it.toUpperCase() }
filtered = samples.findAll { it.contains("1") }
samples.each { println it }

// Maps
config = [input: "data/", threads: 8]
input = config.input               // Dot notation
input = config["input"]            // Bracket notation

// Closures
double = { it * 2 }
result = [1,2,3].collect { it * 2 }  // [2,4,6]

// Control flow
if (x > 5) {
    println "Big"
} else {
    println "Small"
}

for (item in list) {
    println item
}

// You're ready for Nextflow!
```

---

## ✅ Day 2 Completion Checklist

Before marking Day 2 complete, ensure you can:

- [ ] Write string interpolation with `${}`
- [ ] Use `collect` to transform lists
- [ ] Use `findAll` to filter lists
- [ ] Access map values with dot notation
- [ ] Understand the `it` keyword
- [ ] Read a simple Nextflow process
- [ ] Know when to use triple quotes
- [ ] Translate simple Python to Groovy

**Completed Day 2?** Update your `PROGRESS.md`! Tomorrow you'll write real Nextflow code! 🎉

**Your progress**: 2/28 days (7.1%) complete

---

## 🔍 Bonus: Common Gotchas

Watch out for these when starting:

1. **Quotes Matter**:
   ```groovy
   'No interpolation: ${var}'   // Literal string
   "Yes interpolation: ${var}"  // Variable replaced
   ```

2. **Braces Are Your Friend**:
   ```groovy
   "${id}_file.txt"   // ✅ Correct
   "$id_file.txt"     // ❌ Looks for variable "id_file"
   ```

3. **List Methods Are Different**:
   ```groovy
   len(list)       // ❌ Python way
   list.size()     // ✅ Groovy way
   ```

4. **True/False Are Lowercase**:
   ```groovy
   True   // ❌ Python
   true   // ✅ Groovy
   ```

5. **Dot vs Bracket for Maps**:
   ```groovy
   map.simple_key     // ✅ Works
   map."key-with-dash"  // ✅ Dot with quotes for special chars
   map["any_key"]     // ✅ Always works
   ```

---

*Tomorrow: Day 3 - Your First Nextflow Process*

**You're making great progress! See you tomorrow! 🚀**