# Day 2: Groovy Essentials for Nextflow

**Learning Time**: 30 minutes  
**Prerequisites**: Completed Day 1  
**Goal**: Learn just enough Groovy to read and write Nextflow workflows fluently

---

## 📖 Introduction (2 minutes)

Welcome to Day 2! Yesterday you learned *what* Nextflow is and *why* it matters. Today you'll learn *how* to read and write it.

Here's the good news: **You don't need to become a Groovy expert.** Nextflow workflows use only a small subset of Groovy, and you'll learn all of it today.

Think of it like learning to read music. You don't need to know everything about music theory—just enough to read the notes and play the songs you want.

### Why Groovy?

Nextflow is written in Groovy, a dynamic language for the Java Virtual Machine (JVM). But don't worry about Java—Nextflow abstracts away the complex parts. You'll only see the simple, readable parts.

**Key insight:** Groovy syntax looks like a mix of Java and Python, which makes it ideal for bioinformaticians already familiar with Python.

---

## 🎯 Learning Objectives

By the end of Day 2, you should be able to:

1. **Read Groovy code** - Understand any Groovy snippet in a Nextflow workflow
2. **Understand string interpolation** - The most important Groovy feature for Nextflow
3. **Work with data structures** - Lists, maps, and tuples
4. **Recognize closures** - Understand functional programming patterns
5. **Write simple Groovy** - Create basic variables, functions, and expressions
6. **Spot patterns** - Recognize common Groovy idioms in workflows
7. **Debug syntax errors** - Know what went wrong and why

---

## 📚 Key Concepts (20 minutes)

### 1. Groovy vs Python: A Comparison

Both are dynamic languages, but with different syntax. Here's what you need to know:

| Feature | Python | Groovy | Note |
|---------|--------|--------|------|
| **Variables** | `name = "Alice"` | `name = "Alice"` | Same! |
| **Strings** | `f"Hello {name}"` | `"Hello ${name}"` | Use `${}` in Groovy |
| **Lists** | `[1, 2, 3]` | `[1, 2, 3]` | Same syntax |
| **Dicts** | `{"a": 1}` | `[a: 1]` | Groovy uses `[key: value]` |
| **Functions** | `def func():` | `def func() {}` | Similar |
| **Loops** | `for x in items:` | `items.each { x -> ... }` | Different style |
| **Lambda** | `lambda x: x*2` | `{ x -> x*2 }` | Groovy has closures |

**Bottom line:** If you know Python, Groovy will feel familiar but slightly different. The good news is the differences are minimal for what Nextflow uses.

---

### 2. String Interpolation (MOST IMPORTANT!)

String interpolation is **THE** most important Groovy feature for Nextflow. You'll use it constantly.

#### Python Way:
```python
sample_id = "sample_1"
file_path = f"data/{sample_id}.fastq"
print(f"Processing {sample_id} from {file_path}")
```

#### Groovy Way:
```groovy
sample_id = "sample_1"
file_path = "data/${sample_id}.fastq"
println "Processing ${sample_id} from ${file_path}"
```

**Key difference:** Use `${}` instead of `{}` (but inside double quotes)

#### Examples You'll See in Nextflow:

```groovy
// Simple interpolation
sample = "tumor"
command = "bwa mem reference.fa ${sample}.fastq"

// Property access
file = file("data/sample.bam")
message = "File: ${file.baseName} (${file.size()} bytes)"

// Expressions
reads = 1000000
message = "Processed ${reads / 1000000} million reads"

// In script blocks (very common!)
process ALIGN {
    script:
        """
        # Groovy interpolation happens BEFORE running bash
        # So ${sample} is replaced with actual value
        bwa mem ${reference} ${sample}.fastq > ${sample}.bam
        """
}
```

**Critical concept:** Groovy interpolation (`${}`) happens in Nextflow, **before** the script runs. Bash variables need to be escaped with `\${}`.

---

### 3. Variables and Basic Types

Groovy is dynamically typed, like Python. You don't declare types:

```groovy
// Strings
name = "Alice"
path = 'data/sample.fastq'  // Single quotes work too

// Numbers
count = 42                    // Integer
percent = 75.5               // Double
memory = 16                  // Integer

// Booleans
is_valid = true
has_reads = false

// null
result = null              // Like Python's None

// No type declaration needed!
// Groovy figures it out automatically
```

**Python comparison:**
```python
name = "Alice"              # Python
name = "Alice"              # Groovy - same!
```

---

### 4. Lists (Arrays)

Lists are like Python lists:

```groovy
// Creating lists
samples = [1, 2, 3, 4, 5]
tools = ["fastqc", "bwa", "samtools"]
mixed = [1, "text", true, null]  // Can mix types

// Accessing
first = samples[0]                // 1
last = samples[-1]               // 5

// Common operations
samples.size()                    // 5
samples.isEmpty()                 // false
samples.contains(3)               // true

// List methods
doubled = samples.collect { it * 2 }         // [2, 4, 6, 8, 10]
evens = samples.findAll { it % 2 == 0 }     // [2, 4]
sum = samples.sum()                          // 15

// Iteration
samples.each { sample ->
    println "Sample: ${sample}"
}

// Iteration with index
samples.eachWithIndex { sample, index ->
    println "Sample ${index}: ${sample}"
}
```

**Python equivalent:**
```python
samples = [1, 2, 3, 4, 5]
first = samples[0]
doubled = [x * 2 for x in samples]  # Python list comprehension
evens = [x for x in samples if x % 2 == 0]
for sample in samples:
    print(f"Sample: {sample}")
```

---

### 5. Maps (Dictionaries)

Maps are like Python dicts, but with cleaner syntax:

```groovy
// Creating maps
config = [name: "Alice", age: 30, active: true]

// Python would be: config = {"name": "Alice", "age": 30}
// Groovy is cleaner!

// Accessing values
name = config["name"]              // "Alice"
name = config.name                 // "Alice" - cleaner!
age = config.age                   // 30

// Adding/modifying
config.city = "New York"           // Add new key
config["country"] = "USA"          // Alternative syntax

// Checking keys
config.containsKey("name")         // true
"age" in config                    // true

// Getting all keys/values
config.keySet()                    // [name, age, active, city, country]
config.values()                    // [Alice, 30, true, New York, USA]

// Iteration
config.each { key, value ->
    println "${key}: ${value}"
}
```

**Python equivalent:**
```python
config = {"name": "Alice", "age": 30, "active": True}
name = config["name"]              # "Alice"
age = config["age"]                # 30
config["city"] = "New York"
```

---

### 6. Closures (Functions as Values)

Closures are Groovy's version of lambda functions. They're blocks of code you can pass around:

#### Basic Closure Syntax:

```groovy
// A closure is code in curly braces
{ println "Hello!" }

// With parameters
{ name -> println "Hello, ${name}!" }

// With multiple parameters
{ x, y -> x + y }

// The special 'it' parameter (implicit)
{ it * 2 }          // 'it' is the default parameter name
```

#### Using Closures:

```groovy
// Define a closure
double = { x -> x * 2 }

// Call it
result = double(5)                 // 10

// Common with list operations
numbers = [1, 2, 3, 4, 5]

// Transform with collect (map in Python)
doubled = numbers.collect { it * 2 }    // [2, 4, 6, 8, 10]

// Filter
evens = numbers.findAll { it % 2 == 0 } // [2, 4]

// Any/all checks
hasEven = numbers.any { it % 2 == 0 }   // true
allPositive = numbers.all { it > 0 }    // true

// Simple iteration
numbers.each { num ->
    println "Number: ${num}"
}
```

**Python equivalent:**
```python
numbers = [1, 2, 3, 4, 5]
doubled = [x * 2 for x in numbers]      # List comprehension
doubled = list(map(lambda x: x * 2, numbers))  # Using lambda
evens = [x for x in numbers if x % 2 == 0]
evens = list(filter(lambda x: x % 2 == 0, numbers))
```

#### Closures in Nextflow Context:

You'll see closures everywhere in Nextflow channel operations:

```groovy
// Transform files into names
files = Channel.fromPath("data/*.fastq")
names = files.map { file -> file.baseName }

// Filter by file size
large_files = files.filter { file -> file.size() > 1000000 }

// Combine into tuples
paired = files
    .map { file -> [file.baseName, file] }
    .filter { name, file -> file.size() > 0 }
```

---

### 7. String Types and GStrings

Groovy has two types of strings:

```groovy
// Single quotes: NO interpolation
single = 'Hello ${name}'    // Literal: "Hello ${name}"

// Double quotes: Interpolation happens!
double = "Hello ${name}"    // If name="Alice": "Hello Alice"

// This is why Nextflow scripts use double quotes for interpolation
process MY_PROCESS {
    script:
        """
        # Groovy interpolation (double quotes)
        echo "Processing ${sample}"
        
        # Bash variable (single quotes or escape)
        echo "Path is \$PATH"
        """
}
```

---

### 8. Groovy-specific Features You'll See

#### Elvis Operator (?:)

Provide a default if value is null:

```groovy
name = params.name ?: "default"  // Use params.name, or "default" if null
```

#### Safe Navigation (?.)

Prevent errors when accessing properties on null:

```groovy
file = file("data/sample.bam")
size = file?.size()              // Returns null if file is null, not error
```

#### Range

Create sequences:

```groovy
numbers = (1..5)                 // [1, 2, 3, 4, 5]
letters = ('a'..'z')             // ['a', 'b', ..., 'z']

for (i in 1..3) {
    println i
}
```

---

### 9. Methods and Simple Functions

```groovy
// Define a function
def greet(name) {
    return "Hello, ${name}!"
}

// Call it
message = greet("Alice")

// Functions can be one-liners
def double(x) { return x * 2 }

// or without return keyword (implicit return)
def triple(x) { x * 3 }

// With default parameters
def greet(name, greeting = "Hello") {
    return "${greeting}, ${name}!"
}

greet("Alice")                   // "Hello, Alice!"
greet("Bob", "Hi")               // "Hi, Bob!"
```

---

### 10. Control Flow

Same as Python, just different syntax:

```groovy
// If statements
if (x > 10) {
    println "Big"
} else if (x > 5) {
    println "Medium"
} else {
    println "Small"
}

// For loops
for (i in 1..5) {
    println i
}

// While loops
count = 0
while (count < 5) {
    println count
    count++
}

// Switch
switch (value) {
    case "fastq":
        tool = "fastqc"
        break
    case "bam":
        tool = "samtools"
        break
    default:
        tool = "unknown"
}

// Try-catch
try {
    result = doSomething()
} catch (Exception e) {
    println "Error: ${e.message}"
}
```

---

## 🔗 Python to Groovy Cheat Sheet

| Python | Groovy | Notes |
|--------|--------|-------|
| `name = "Alice"` | `name = "Alice"` | Same |
| `f"Hello {name}"` | `"Hello ${name}"` | Use `${}` for interpolation |
| `if x > 5:` | `if (x > 5) {` | Parentheses required |
| `for x in items:` | `items.each { x -> ... }` | Groovy uses closures |
| `[1,2,3]` | `[1, 2, 3]` | Same |
| `{"a": 1}` | `[a: 1]` | Different syntax |
| `lambda x: x*2` | `{ x -> x*2 }` | Closures vs lambdas |
| `dict.keys()` | `map.keySet()` | Different method names |
| `list[0]` | `list[0]` | Same indexing |
| `len(list)` | `list.size()` | Different method names |
| `print()` | `println` | Groovy-style |
| `None` | `null` | Different keyword |
| `True/False` | `true/false` | Lowercase |
| `def func():` | `def func() {}` | Similar |

---

## 💻 Real Nextflow Examples

### Example 1: Understanding a Simple Process

```groovy
process QUALITY_CHECK {
    // Input declaration
    input:
        path fastq_file          // Groovy: Variable with type
        val sample_id            // Groovy: Another variable
    
    // Output declaration
    output:
        path "*.html"
    
    // Script: The actual work
    script:
        """
        # This is bash, but Groovy variables are interpolated first!
        fastqc --threads 4 ${fastq_file}
        """
}
```

**What's happening:**
- `fastq_file` and `sample_id` are Groovy variables
- Inside the `script:` block with triple quotes, `${fastq_file}` gets replaced with the actual value
- The `*.html` is a glob pattern (string)

### Example 2: Using Parameters

```groovy
// Define parameters with defaults
params {
    input_dir = "data/"
    output_dir = "results/"
    quality_threshold = 30
}

// Use them in a process
process FILTER_READS {
    input:
        path reads
    output:
        path "*.filtered.fastq"
    script:
        """
        # Groovy interpolation happens
        seqtk seq -q ${params.quality_threshold} ${reads} > output.filtered.fastq
        """
}
```

**Key points:**
- `params.quality_threshold` - Accessing a map value with dot notation
- `${params.quality_threshold}` - Interpolating into a string
- `params` is a special Nextflow map

### Example 3: Lists and Channels

```groovy
// Python-style list
tools = ["fastqc", "bwa", "samtools"]

// Using with closure
tools.each { tool ->
    println "Tool: ${tool}"
}

// Transform with closure
uppercase_tools = tools.collect { tool -> tool.toUpperCase() }
// Result: ["FASTQC", "BWA", "SAMTOOLS"]

// Channel with list
Channel.from(tools)
    .map { tool -> "${tool}_v1.0" }  // Add version
    .view()                           // Print each item
```

### Example 4: Maps in Workflows

```groovy
// Define sample metadata
samples = [
    [id: "sample1", fastq: "sample1.fastq", condition: "control"],
    [id: "sample2", fastq: "sample2.fastq", condition: "treatment"],
]

// Use with closure
samples.each { sample ->
    println "Processing ${sample.id} (${sample.condition})"
}

// Transform
sample_ids = samples.collect { sample -> sample.id }
// Result: ["sample1", "sample2"]

// Filter
controls = samples.findAll { sample -> sample.condition == "control" }
```

### Example 5: Conditional Logic

```groovy
process ALIGN {
    input:
        val sample_id
        path fastq
        path reference
    
    script:
        // Groovy conditional - chooses script at runtime
        if (params.aligner == "bwa") {
            """
            bwa mem ${reference} ${fastq} > ${sample_id}.sam
            """
        } else if (params.aligner == "bowtie2") {
            """
            bowtie2 -x ${reference} -U ${fastq} -S ${sample_id}.sam
            """
        } else {
            error "Unknown aligner: ${params.aligner}"
        }
}
```

**Key insight:** The `script:` block is evaluated as Groovy at *workflow definition time*, not at *execution time*. So the `if` chooses which script to run.

---

## 🧠 Common Groovy Patterns in Nextflow

### Pattern 1: Safe Access

```groovy
// Safe navigation - returns null if not found
size = file?.size()

// Elvis operator - provide default
name = params.name ?: "default"
```

### Pattern 2: Multi-value Returns

```groovy
// Return tuple from closure
data = files.map { file ->
    [file.baseName, file]  // Return list as tuple
}
```

### Pattern 3: Method Chaining

```groovy
// Chain operations on channels/lists
result = samples
    .map { s -> [s.id, s.fastq] }           // Transform
    .filter { id, fastq -> fastq.size() > 0 }  // Filter
    .collect()                               // Gather results
```

### Pattern 4: Closure with Multiple Parameters

```groovy
// Two parameters
samples.map { sample, index ->
    "${index}: ${sample}"
}

// Destructuring in closure
paired_files.each { (read1, read2) ->
    println "Processing pair: ${read1} and ${read2}"
}
```

---

## 📊 Groovy Cheat Sheet Reference

### Variables
```groovy
name = "Alice"          // String
age = 30                // Integer
price = 19.99           // Double
active = true           // Boolean
nothing = null          // null
```

### Strings
```groovy
double_quotes = "Hello ${name}"      // Interpolation
single_quotes = 'Hello ${name}'      // Literal
multiline = """
Multiple
lines
here
"""
```

### Collections
```groovy
list = [1, 2, 3]                    // List
map = [a: 1, b: 2]                  // Map
tuple = [1, "text"]                 // List as tuple
range = (1..5)                      // Range [1,2,3,4,5]
```

### Accessing
```groovy
list[0]                             // First item
map.key                             // Map value
map["key"]                          // Alternative
list.last()                         // Last item
```

### Operations
```groovy
"text".toUpperCase()                // UPPERCASE
"text".size()                       // 4
[1,2,3].size()                      // 3
[1,2,3].sum()                       // 6
file.baseName                       // Filename without extension
```

### Closures
```groovy
{ x -> x * 2 }                      // Simple closure
{ it * 2 }                          // Using 'it'
items.collect { item -> ... }       // Transform
items.findAll { item -> ... }       // Filter
items.each { item -> ... }          // Iterate
```

### Control Flow
```groovy
if (x > 5) { ... }                  // If statement
for (i in 1..5) { ... }             // Loop
while (true) { ... }                // While loop
try { ... } catch (e) { ... }       // Error handling
```

---

## 🚨 Common Groovy Gotchas for Python Developers

### Gotcha 1: Dollar Sign Interpolation

**Wrong:**
```groovy
message = f"Hello {name}"           // Python syntax - doesn't work in Groovy!
message = "Hello {name}"            // No interpolation
```

**Right:**
```groovy
message = "Hello ${name}"           // Correct Groovy syntax
```

### Gotcha 2: Map Syntax

**Wrong:**
```groovy
config = {"name": "Alice"}          // Python dict syntax
```

**Right:**
```groovy
config = [name: "Alice"]            // Groovy map syntax
```

### Gotcha 3: Loops

**Wrong:**
```groovy
for item in items { }               // Python syntax
```

**Right:**
```groovy
for (item in items) { }             // Groovy syntax with parentheses
items.each { item -> }              // Groovy-style with closure
```

### Gotcha 4: Null vs Uninitialized

```groovy
x = null                            // Explicitly null
// x is not defined - will error if accessed

// Use Elvis operator for defaults
y = params.value ?: "default"
```

### Gotcha 5: String vs Bash Variables in Scripts

**In Nextflow:**
```groovy
process MY_PROCESS {
    input:
        val sample_id
    
    script:
        """
        # Groovy interpolation (processed by Nextflow first)
        SAMPLE=${sample_id}
        
        # Bash variable (processed by bash, \$ escapes the $)
        echo "Path: \$PATH"
        
        # Bash variable expansion (after substitution)
        echo "Sample: \$SAMPLE"
        """
}
```

When `sample_id = "sample1"`:
```bash
# What actually runs:
SAMPLE=sample1
echo "Path: $PATH"
echo "Sample: $SAMPLE"
```

---

## ✅ Completion Checklist

Review these items to confirm you understand:

- [ ] I understand string interpolation with `${variable}`
- [ ] I can read and write Groovy lists
- [ ] I can read and write Groovy maps
- [ ] I understand closures with `{ ... }`
- [ ] I know when to use single vs double quotes
- [ ] I recognize common Groovy patterns
- [ ] I can spot Groovy errors and guess why they happened
- [ ] I feel comfortable reading Nextflow workflows
- [ ] I could write a simple Groovy script

---

## 🔑 Key Takeaways

### The Essential Groovy for Nextflow

1. **String Interpolation** - `"Hello ${name}"` - You'll use this constantly
2. **Lists** - `[1, 2, 3]` - Same as Python
3. **Maps** - `[a: 1, b: 2]` - Like Python dicts
4. **Closures** - `{ x -> x * 2 }` - Transform and filter data
5. **Safe Access** - `obj?.property` - Don't crash on null
6. **Method Chaining** - `.map().filter().collect()` - Compose operations

### What You DON'T Need to Know

- ❌ Advanced Groovy features
- ❌ Java integration
- ❌ Metaclasses or metaprogramming
- ❌ Groovy-specific libraries
- ❌ Complex OOP concepts

You only need the small subset used in Nextflow workflows!

---

## 🎯 How This Prepares You for Day 3

Tomorrow (Day 3), you'll write your first Nextflow process. You'll use:
- String interpolation (to build commands)
- Maps (to organize inputs/outputs)
- Closures (to handle data)
- Basic control flow (if/for statements)

Everything you learned today will make Day 3 straightforward.

---

## 📚 Quick Reference: Python → Groovy

When you see Nextflow code and think "what is that?", check this table:

| Pattern | Example | Explanation |
|---------|---------|-------------|
| String with value | `"Hello ${name}"` | Interpolation (use `${}`) |
| List | `[1, 2, 3]` | Literal list |
| Map | `[key: value]` | Literal map |
| Transform list | `items.collect { it * 2 }` | Apply closure to each |
| Filter list | `items.findAll { it > 5 }` | Keep matching items |
| Loop | `items.each { item -> }` | Execute for each |
| Access map | `map.key` or `map["key"]` | Get value |
| Default if null | `value ?: "default"` | Elvis operator |
| Safe access | `obj?.property` | Returns null, not error |

---

## 🎓 Self-Assessment

**Can you do these?**

1. **Read code:** Look at a Nextflow process and explain what the Groovy parts do
2. **Understand interpolation:** Explain what `"File: ${file.baseName}.bam"` produces
3. **Write closures:** Write a closure that transforms a list
4. **Use maps:** Create a map and access its values
5. **Identify errors:** Look at wrong Groovy syntax and spot the mistake

If you can do all 5, you're ready for Day 3!

---

## 🚀 Ready for Exercises?

Turn the page to try the hands-on exercises. Don't worry if you make mistakes—that's how you learn!

The exercises will have you:
1. Reading Groovy code
2. Predicting output
3. Writing simple Groovy
4. Debugging syntax errors
5. Understanding closure patterns

Remember: The goal isn't perfection. It's familiarity and confidence. By the end, you'll feel comfortable with Groovy in a Nextflow context.

---

*This is Day 2 of 28. You're building the foundation for writing production-ready Nextflow workflows!*
