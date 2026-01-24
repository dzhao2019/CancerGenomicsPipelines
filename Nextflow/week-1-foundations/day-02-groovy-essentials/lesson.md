\# Day 2: Groovy Essentials for Nextflow



\*\*Learning Time\*\*: 30 minutes  

\*\*Prerequisites\*\*: Day 1 completed, Basic Python knowledge  

\*\*Goal\*\*: Read and write Groovy code comfortably for Nextflow workflows



---



\## 📖 Introduction (3 minutes)



Welcome to Day 2! Today you'll learn Groovy—but don't worry, you're not becoming a Groovy expert. You're learning just enough to read and write Nextflow workflows comfortably.



\*\*The Good News\*\*: If you know Python, you already understand 80% of what you need. Groovy is like Python's Java-based cousin—similar ideas, slightly different syntax.



\### What You'll Learn Today



\- Core Groovy syntax you'll use in Nextflow

\- String interpolation (the `${}` pattern you'll see everywhere)

\- Closures (Groovy's version of lambda functions)

\- Lists and Maps (similar to Python lists and dicts)

\- Key differences from Python that matter

\- What you DON'T need to learn (most of Groovy!)



\### What You WON'T Learn



You're not learning:

\- ❌ Object-oriented Groovy programming

\- ❌ Advanced Groovy features

\- ❌ Groovy frameworks

\- ❌ Java interoperability details



You're learning:

\- ✅ Reading Nextflow workflow syntax

\- ✅ Writing simple Groovy expressions

\- ✅ Working with data structures in Nextflow

\- ✅ Understanding patterns you'll see daily



---



\## 🎯 The Big Picture: Why Groovy? (2 minutes)



\*\*Why not just use Python for Nextflow?\*\*



Nextflow chose Groovy because:

1\. \*\*JVM-based\*\*: Runs on the Java Virtual Machine (portable, robust)

2\. \*\*Dynamic\*\*: Similar feel to Python (no type declarations needed)

3\. \*\*Concise\*\*: Less verbose than Java

4\. \*\*DSL-friendly\*\*: Great for creating Domain Specific Languages (which Nextflow is)



\*\*For you\*\*: It means learning some new syntax, but the logic is familiar.



---



\## 🔑 Key Concepts with Examples (15 minutes)



\### 1. Variables and Basic Types (2 minutes)



\*\*Python\*\*:

```python

name = "Alice"

age = 30

height = 5.6

is\_scientist = True

```



\*\*Groovy\*\*:

```groovy

name = "Alice"

age = 30

height = 5.6

is\_scientist = true  // lowercase true/false

```



\*\*Key Differences\*\*:

\- Groovy uses `true/false` (lowercase) instead of `True/False`

\- Groovy uses `null` instead of `None`

\- Semicolons are optional (usually omitted)

\- Type declarations are optional



```groovy

// All of these work in Groovy

def name = "Alice"        // 'def' means dynamic type

String name = "Alice"     // Explicit type (optional)

name = "Alice"            // Type inferred (most common in Nextflow)

```



\*\*In Nextflow, you'll see\*\*:

```groovy

params.input = "data/"

params.output = "results/"

params.threads = 8

```



\### 2. String Interpolation - THE MOST IMPORTANT PATTERN (3 minutes)



This is crucial because you'll use it constantly in Nextflow!



\*\*Python\*\*:

```python

name = "Alice"

age = 30



\# F-strings (modern Python)

message = f"Hello, {name}! You are {age} years old."



\# Format method

message = "Hello, {}! You are {} years old.".format(name, age)



\# Old style

message = "Hello, %s! You are %d years old." % (name, age)

```



\*\*Groovy\*\*:

```groovy

name = "Alice"

age = 30



// String interpolation with ${}

message = "Hello, ${name}! You are ${age} years old."



// For simple variables, you can skip the braces

message = "Hello, $name! You are $age years old."

```



\*\*When to use braces\*\*:

```groovy

// ALWAYS use braces for:

// 1. Expressions

output = "Result: ${x + y}"



// 2. Property access

output = "File: ${file.name}"



// 3. When variable name is followed by more text

output = "sample\_${id}\_processed.bam"  // Correct

output = "sample\_$id\_processed.bam"    // Wrong! Looks for variable "id\_processed"

```



\*\*In Nextflow scripts\*\* (you'll see this everywhere):

```groovy

process align {

&nbsp;   input:

&nbsp;   val sample\_id

&nbsp;   path reads

&nbsp;   

&nbsp;   output:

&nbsp;   path "${sample\_id}.bam"

&nbsp;   

&nbsp;   script:

&nbsp;   """

&nbsp;   bwa mem reference.fa ${reads} > ${sample\_id}.bam

&nbsp;   samtools sort ${sample\_id}.bam -o ${sample\_id}.sorted.bam

&nbsp;   """

}

```



\*\*Pro Tip\*\*: In Nextflow `script` blocks (triple quotes), you use `${variable}` to reference Groovy variables. This is THE pattern you'll use most!



\### 3. Lists (Like Python Lists) (2 minutes)



\*\*Python\*\*:

```python

\# Creating lists

samples = \["sample1", "sample2", "sample3"]

numbers = \[1, 2, 3, 4, 5]



\# Accessing elements

first = samples\[0]

last = samples\[-1]



\# List operations

samples.append("sample4")

length = len(samples)



\# List comprehension

doubled = \[x \* 2 for x in numbers]

```



\*\*Groovy\*\*:

```groovy

// Creating lists

samples = \["sample1", "sample2", "sample3"]

numbers = \[1, 2, 3, 4, 5]



// Accessing elements

first = samples\[0]

last = samples\[-1]  // Yes, negative indexing works!



// List operations

samples << "sample4"  // << is like append

samples.add("sample5")  // .add() also works

length = samples.size()



// Using collect (like map)

doubled = numbers.collect { it \* 2 }

```



\*\*Key Differences\*\*:

\- `<<` operator adds to list (common in Groovy)

\- `size()` instead of `len()`

\- `collect {}` is like Python's list comprehension



\*\*Common Groovy List Methods\*\*:

```groovy

samples = \["sample1", "sample2", "sample3"]



// Filter

long\_names = samples.findAll { it.length() > 7 }



// Map/transform

uppercase = samples.collect { it.toUpperCase() }



// Each (like for loop)

samples.each { println it }



// Join

text = samples.join(", ")  // "sample1, sample2, sample3"

```



\*\*In Nextflow\*\*:

```groovy

// You'll often see

Channel.from(\["sample1", "sample2", "sample3"])



// Or using ranges

Channel.from(1..10)  // Numbers 1 through 10

```



\### 4. Maps (Like Python Dictionaries) (2 minutes)



\*\*Python\*\*:

```python

\# Creating dictionaries

config = {

&nbsp;   "input": "data/",

&nbsp;   "output": "results/",

&nbsp;   "threads": 8

}



\# Accessing values

input\_dir = config\["input"]

threads = config.get("threads", 4)



\# Adding/updating

config\["genome"] = "hg38"

```



\*\*Groovy\*\*:

```groovy

// Creating maps

config = \[

&nbsp;   input: "data/",

&nbsp;   output: "results/",

&nbsp;   threads: 8

]



// Accessing values (multiple ways)

input\_dir = config\["input"]     // Bracket notation

input\_dir = config.input        // Dot notation

threads = config.get("threads", 4)



// Adding/updating

config.genome = "hg38"

config\["reference"] = "ref.fa"

```



\*\*Key Differences\*\*:

\- No quotes needed for keys (unless they have spaces/special chars)

\- Can use dot notation for access

\- Square brackets `\[]` for maps, not curly braces `{}`



\*\*In Nextflow Parameters\*\*:

```groovy

// Very common pattern

params.input = "data/\*.fastq"

params.output = "results/"

params.genome = "hg38"

params.threads = 8



// Accessing in processes

process align {

&nbsp;   cpus params.threads

&nbsp;   

&nbsp;   script:

&nbsp;   """

&nbsp;   bwa mem -t ${params.threads} reference.fa reads.fastq

&nbsp;   """

}

```



\### 5. Closures - Groovy's Lambda Functions (3 minutes)



This is where Groovy differs most from Python, but it's actually quite elegant!



\*\*Python Lambda\*\*:

```python

\# Simple lambda

double = lambda x: x \* 2



\# With map

numbers = \[1, 2, 3, 4, 5]

doubled = list(map(lambda x: x \* 2, numbers))



\# With filter

evens = list(filter(lambda x: x % 2 == 0, numbers))

```



\*\*Groovy Closures\*\*:

```groovy

// Closure (like lambda)

double = { x -> x \* 2 }



// With collect (like map)

numbers = \[1, 2, 3, 4, 5]

doubled = numbers.collect { x -> x \* 2 }



// Implicit 'it' parameter (super common!)

doubled = numbers.collect { it \* 2 }



// With findAll (like filter)

evens = numbers.findAll { it % 2 == 0 }

```



\*\*The `it` keyword\*\*: When a closure has one parameter, you can use `it` instead of naming it.



```groovy

// These are equivalent:

samples.collect { sample -> sample.toUpperCase() }

samples.collect { it.toUpperCase() }



// These are equivalent:

numbers.findAll { n -> n > 5 }

numbers.findAll { it > 5 }

```



\*\*Multi-line Closures\*\*:

```groovy

// For complex operations

results = samples.collect { sample ->

&nbsp;   def trimmed = sample.trim()

&nbsp;   def upper = trimmed.toUpperCase()

&nbsp;   return "${upper}\_processed"

}



// Last expression is automatically returned

results = samples.collect { sample ->

&nbsp;   def trimmed = sample.trim()

&nbsp;   trimmed.toUpperCase()  // Returned automatically

}

```



\*\*In Nextflow\*\* (you'll see closures everywhere):

```groovy

// Transforming channel data

Channel

&nbsp;   .fromPath("\*.fastq")

&nbsp;   .map { file -> \[file.baseName, file] }

&nbsp;   .view { name, file -> "Processing: ${name}" }



// Using 'it' for simple cases

Channel

&nbsp;   .fromPath("\*.fastq")

&nbsp;   .map { it.baseName }

&nbsp;   .view { "Sample: ${it}" }

```



\### 6. Strings: Single vs Double Quotes (2 minutes)



\*\*Python\*\*: Single and double quotes are mostly interchangeable

```python

name = 'Alice'  # Same as "Alice"

message = f"Hello, {name}"  # F-string for interpolation

```



\*\*Groovy\*\*: Quote type matters!

```groovy

// Single quotes: NO interpolation (literal string)

message = 'Hello, ${name}'  // Prints: Hello, ${name}



// Double quotes: YES interpolation

message = "Hello, ${name}"  // Prints: Hello, Alice



// Triple quotes: Multi-line with interpolation

script = """

&nbsp;   echo "Processing ${sample}"

&nbsp;   bwa mem ${reference} ${reads}

&nbsp;   """

```



\*\*In Nextflow\*\*, you'll see:

```groovy

process example {

&nbsp;   script:

&nbsp;   """

&nbsp;   # This is a multi-line shell script

&nbsp;   # ${variables} are interpolated from Groovy

&nbsp;   echo "Sample: ${sample\_id}"

&nbsp;   fastqc ${input\_file}

&nbsp;   """

}

```



\*\*Important Rule\*\*: 

\- `"""triple quotes"""` for multi-line scripts (Nextflow script blocks)

\- `"double quotes"` for single-line strings with interpolation

\- `'single quotes'` for literal strings (rare in Nextflow)



\### 7. Control Flow - Familiar Territory (2 minutes)



\*\*If/Else\*\*:



Python:

```python

if score > 90:

&nbsp;   grade = "A"

elif score > 80:

&nbsp;   grade = "B"

else:

&nbsp;   grade = "C"

```



Groovy:

```groovy

if (score > 90) {

&nbsp;   grade = "A"

} else if (score > 80) {

&nbsp;   grade = "B"

} else {

&nbsp;   grade = "C"

}

```



\*\*For Loops\*\*:



Python:

```python

for sample in samples:

&nbsp;   print(sample)



for i in range(10):

&nbsp;   print(i)

```



Groovy:

```groovy

// For-each style

for (sample in samples) {

&nbsp;   println sample

}



// Using each (more common in Groovy)

samples.each { sample ->

&nbsp;   println sample

}



// Ranges

for (i in 0..9) {

&nbsp;   println i

}



(0..9).each { println it }

```



\*\*In Nextflow\*\*, you typically use channels and operators instead of explicit loops, but you'll see these patterns in script blocks:



```groovy

process example {

&nbsp;   script:

&nbsp;   """

&nbsp;   # Shell script can have normal bash loops

&nbsp;   for file in \*.bam; do

&nbsp;       echo \\$file  # Note: \\$ for bash variables

&nbsp;   done

&nbsp;   

&nbsp;   # But Groovy ${} for Nextflow variables

&nbsp;   echo "Sample ID: ${sample\_id}"

&nbsp;   """

}

```



\### 8. What You DON'T Need to Know (1 minute)



\*\*Relax—you can skip\*\*:

\- Classes and object-oriented programming

\- Groovy's advanced metaprogramming

\- Gradle and Groovy build tools

\- Most of the Groovy standard library

\- Type systems and generics

\- Exception handling details



\*\*Focus on\*\*:

\- String interpolation `${}`

\- Lists and Maps

\- Closures with `it`

\- Reading Nextflow examples



---



\## 🔗 Python Connection: Side-by-Side Comparison (3 minutes)



Let's see a complete example in both languages:



\### Task: Process a list of samples



\*\*Python Version\*\*:

```python

samples = \["sample1", "sample2", "sample3"]



\# Transform: add suffix

processed = \[f"{s}\_processed.bam" for s in samples]



\# Filter: only long names

long\_names = \[s for s in samples if len(s) > 7]



\# Print each

for sample in samples:

&nbsp;   print(f"Processing {sample}")



\# Dictionary of sample info

sample\_info = {

&nbsp;   "sample1": {"type": "tumor", "batch": 1},

&nbsp;   "sample2": {"type": "normal", "batch": 1},

&nbsp;   "sample3": {"type": "tumor", "batch": 2}

}



\# Access nested data

for sample, info in sample\_info.items():

&nbsp;   print(f"{sample}: {info\['type']} from batch {info\['batch']}")

```



\*\*Groovy Version\*\*:

```groovy

samples = \["sample1", "sample2", "sample3"]



// Transform: add suffix

processed = samples.collect { "${it}\_processed.bam" }



// Filter: only long names

long\_names = samples.findAll { it.length() > 7 }



// Print each

samples.each { sample ->

&nbsp;   println "Processing ${sample}"

}



// Map of sample info

sample\_info = \[

&nbsp;   sample1: \[type: "tumor", batch: 1],

&nbsp;   sample2: \[type: "normal", batch: 1],

&nbsp;   sample3: \[type: "tumor", batch: 2]

]



// Access nested data

sample\_info.each { sample, info ->

&nbsp;   println "${sample}: ${info.type} from batch ${info.batch}"

}

```



\### Common Patterns Translation Table



| Task | Python | Groovy |

|------|--------|--------|

| String interpolation | `f"{var}"` | `"${var}"` or `"$var"` |

| List comprehension | `\[x\*2 for x in list]` | `list.collect { it\*2 }` |

| Filter list | `\[x for x in list if x>5]` | `list.findAll { it>5 }` |

| For each | `for x in list: print(x)` | `list.each { println it }` |

| Dictionary | `{"key": "value"}` | `\[key: "value"]` |

| Access dict | `dict\["key"]` or `dict.get("key")` | `dict\["key"]` or `dict.key` |

| Length | `len(list)` | `list.size()` |

| Append to list | `list.append(x)` | `list << x` or `list.add(x)` |

| Join strings | `", ".join(list)` | `list.join(", ")` |

| True/False | `True/False` | `true/false` |

| None/null | `None` | `null` |



---



\## 💻 Hands-On Exercises (10 minutes)



\### Exercise 1: String Interpolation Practice (3 minutes)



Convert these Python strings to Groovy with proper interpolation:



\*\*Python\*\*:

```python

sample\_id = "patient001"

replicate = 2

file\_type = "fastq"



\# Task 1: Create filename

filename = f"{sample\_id}\_rep{replicate}.{file\_type}"



\# Task 2: Create command

threads = 8

command = f"fastqc -t {threads} {filename}"



\# Task 3: Create output path

output\_dir = "/results"

output\_path = f"{output\_dir}/{sample\_id}/qc\_rep{replicate}.html"

```



\*\*Your Groovy Version\*\*:

```groovy

sample\_id = "patient001"

replicate = 2

file\_type = "fastq"



// Task 1: Create filename



// Task 2: Create command



// Task 3: Create output path

```



<details>

<summary>Click to see solution</summary>



```groovy

sample\_id = "patient001"

replicate = 2

file\_type = "fastq"



// Task 1: Create filename

filename = "${sample\_id}\_rep${replicate}.${file\_type}"

// Result: "patient001\_rep2.fastq"



// Task 2: Create command

threads = 8

command = "fastqc -t ${threads} ${filename}"

// Result: "fastqc -t 8 patient001\_rep2.fastq"



// Task 3: Create output path

output\_dir = "/results"

output\_path = "${output\_dir}/${sample\_id}/qc\_rep${replicate}.html"

// Result: "/results/patient001/qc\_rep2.html"



// Alternative for simple cases (no braces needed):

filename = "${sample\_id}\_rep${replicate}.$file\_type"  // Works!

// But safer to always use braces for consistency

```



\*\*Key Points\*\*:

\- Use `${}` for all variable interpolation

\- Braces are required when variable name runs into other text

\- Double quotes enable interpolation, single quotes don't

</details>



\### Exercise 2: Lists and Closures (4 minutes)



Convert these Python operations to Groovy:



\*\*Python\*\*:

```python

samples = \["sample1", "sample2", "sample3", "sample4", "sample5"]



\# Task 1: Create list of filenames like "sample1.bam"

bam\_files = \[f"{s}.bam" for s in samples]



\# Task 2: Filter samples that contain "1" or "2"

early\_samples = \[s for s in samples if "1" in s or "2" in s]



\# Task 3: Print each sample with index

for i, sample in enumerate(samples):

&nbsp;   print(f"Sample {i}: {sample}")



\# Task 4: Create uppercase versions

upper\_samples = \[s.upper() for s in samples]



\# Task 5: Get samples longer than 7 characters

long\_names = \[s for s in samples if len(s) > 7]

```



\*\*Your Groovy Version\*\*:

```groovy

samples = \["sample1", "sample2", "sample3", "sample4", "sample5"]



// Task 1: Create list of filenames like "sample1.bam"



// Task 2: Filter samples that contain "1" or "2"



// Task 3: Print each sample with index



// Task 4: Create uppercase versions



// Task 5: Get samples longer than 7 characters

```



<details>

<summary>Click to see solution</summary>



```groovy

samples = \["sample1", "sample2", "sample3", "sample4", "sample5"]



// Task 1: Create list of filenames like "sample1.bam"

bam\_files = samples.collect { "${it}.bam" }

// Result: \["sample1.bam", "sample2.bam", "sample3.bam", ...]



// Alternative with explicit parameter name:

bam\_files = samples.collect { sample -> "${sample}.bam" }



// Task 2: Filter samples that contain "1" or "2"

early\_samples = samples.findAll { it.contains("1") || it.contains("2") }

// Result: \["sample1", "sample2"]



// Task 3: Print each sample with index

samples.eachWithIndex { sample, i ->

&nbsp;   println "Sample ${i}: ${sample}"

}

// Note: eachWithIndex gives you both element and index



// Task 4: Create uppercase versions

upper\_samples = samples.collect { it.toUpperCase() }

// Result: \["SAMPLE1", "SAMPLE2", "SAMPLE3", ...]



// Task 5: Get samples longer than 7 characters

long\_names = samples.findAll { it.length() > 7 }

// Result: \[] (none are longer than 7 in this case)



// If we had longer names:

test = \["sample1", "sample2", "verylongsample"]

long\_names = test.findAll { it.length() > 7 }

// Result: \["verylongsample"]

```



\*\*Groovy Closure Patterns\*\*:

\- `collect { }` = transform each element (like Python map)

\- `findAll { }` = filter elements (like Python filter)

\- `each { }` = iterate (like Python for loop)

\- `eachWithIndex { element, index -> }` = iterate with index

\- Use `it` for simple single-parameter closures

\- Use named parameters `{ x -> }` for clarity with multiple operations

</details>



\### Exercise 3: Working with Maps (3 minutes)



Convert this Python dictionary manipulation to Groovy:



\*\*Python\*\*:

```python

\# Sample metadata

samples = {

&nbsp;   "sample1": {

&nbsp;       "type": "tumor",

&nbsp;       "patient": "P001",

&nbsp;       "batch": 1

&nbsp;   },

&nbsp;   "sample2": {

&nbsp;       "type": "normal",

&nbsp;       "patient": "P001",

&nbsp;       "batch": 1

&nbsp;   },

&nbsp;   "sample3": {

&nbsp;       "type": "tumor",

&nbsp;       "patient": "P002",

&nbsp;       "batch": 2

&nbsp;   }

}



\# Task 1: Access sample1's type

sample1\_type = samples\["sample1"]\["type"]



\# Task 2: Add a new sample

samples\["sample4"] = {

&nbsp;   "type": "normal",

&nbsp;   "patient": "P002",

&nbsp;   "batch": 2

}



\# Task 3: Print all tumor samples

for sample, info in samples.items():

&nbsp;   if info\["type"] == "tumor":

&nbsp;       print(f"{sample} is from patient {info\['patient']}")



\# Task 4: Create list of all patients

patients = \[info\["patient"] for sample, info in samples.items()]

unique\_patients = list(set(patients))

```



\*\*Your Groovy Version\*\*:

```groovy

// Sample metadata

samples = \[

&nbsp;   sample1: \[

&nbsp;       type: "tumor",

&nbsp;       patient: "P001",

&nbsp;       batch: 1

&nbsp;   ],

&nbsp;   sample2: \[

&nbsp;       type: "normal",

&nbsp;       patient: "P001",

&nbsp;       batch: 1

&nbsp;   ],

&nbsp;   sample3: \[

&nbsp;       type: "tumor",

&nbsp;       patient: "P002",

&nbsp;       batch: 2

&nbsp;   ]

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

samples = \[

&nbsp;   sample1: \[

&nbsp;       type: "tumor",

&nbsp;       patient: "P001",

&nbsp;       batch: 1

&nbsp;   ],

&nbsp;   sample2: \[

&nbsp;       type: "normal",

&nbsp;       patient: "P001",

&nbsp;       batch: 1

&nbsp;   ],

&nbsp;   sample3: \[

&nbsp;       type: "tumor",

&nbsp;       patient: "P002",

&nbsp;       batch: 2

&nbsp;   ]

]



// Task 1: Access sample1's type (multiple ways!)

sample1\_type = samples.sample1.type           // Dot notation (cleanest)

sample1\_type = samples\["sample1"]\["type"]     // Bracket notation

sample1\_type = samples.sample1\["type"]        // Mixed

// Result: "tumor"



// Task 2: Add a new sample

samples.sample4 = \[

&nbsp;   type: "normal",

&nbsp;   patient: "P002",

&nbsp;   batch: 2

]

// Or using brackets:

samples\["sample4"] = \[type: "normal", patient: "P002", batch: 2]



// Task 3: Print all tumor samples

samples.each { sample, info ->

&nbsp;   if (info.type == "tumor") {

&nbsp;       println "${sample} is from patient ${info.patient}"

&nbsp;   }

}

// Output:

// sample1 is from patient P001

// sample3 is from patient P002



// Task 4: Create list of all patients and get unique ones

patients = samples.collect { sample, info -> info.patient }

// Result: \["P001", "P001", "P002", "P002"]



unique\_patients = patients.unique()

// Result: \["P001", "P002"]



// Or in one line:

unique\_patients = samples.collect { sample, info -> info.patient }.unique()

```



\*\*Key Map Patterns\*\*:

\- Access with dots: `map.key.subkey`

\- Access with brackets: `map\["key"]\["subkey"]`

\- Iterate: `map.each { key, value -> }`

\- Transform values: `map.collect { key, value -> }`

\- No quotes needed for keys (unless special characters)

</details>



---



\## 🤔 Reflection Activity (5 minutes)



\### Question 1: Groovy vs Python Comfort Level



Look back at the exercises. Which Groovy concepts feel:



\*\*Comfortable\*\* (similar to Python):

\- \[ ] Variables and basic types

\- \[ ] Lists

\- \[ ] If/else statements

\- \[ ] For loops



\*\*Slightly Different\*\* (need practice):

\- \[ ] String interpolation with `${}`

\- \[ ] Maps (dictionaries)

\- \[ ] Dot notation for map access

\- \[ ] Method names (size vs len)



\*\*New and Different\*\* (need more practice):

\- \[ ] Closures and `it`

\- \[ ] `collect {}` and `findAll {}`

\- \[ ] The `<<` operator

\- \[ ] Triple-quoted strings



\### Question 2: Pattern Practice



Write down these common patterns in your own words:



\*\*String with variable\*\*:

```groovy

// How would you write: "sample\_{id}\_processed.bam"?

// Your answer:

```



\*\*Transform a list\*\*:

```groovy

// How would you convert \["a", "b", "c"] to \["A", "B", "C"]?

// Your answer:

```



\*\*Filter a list\*\*:

```groovy

// How would you get numbers > 5 from \[1,2,3,6,7,8]?

// Your answer:

```



<details>

<summary>Check your answers</summary>



```groovy

// String with variable

filename = "${sample\_id}\_processed.bam"



// Transform a list

uppercase = \["a", "b", "c"].collect { it.toUpperCase() }



// Filter a list

big\_numbers = \[1,2,3,6,7,8].findAll { it > 5 }

// Result: \[6, 7, 8]

```

</details>



\### Question 3: Real Nextflow Example



Here's a real Nextflow snippet. Can you understand it now?



```groovy

process trimReads {

&nbsp;   input:

&nbsp;   tuple val(sample\_id), path(reads)

&nbsp;   

&nbsp;   output:

&nbsp;   tuple val(sample\_id), path("${sample\_id}\_trimmed.fastq")

&nbsp;   

&nbsp;   script:

&nbsp;   """

&nbsp;   trimmomatic SE ${reads} ${sample\_id}\_trimmed.fastq LEADING:3 TRAILING:3

&nbsp;   """

}



workflow {

&nbsp;   samples = Channel.fromFilePairs("data/\*\_{1,2}.fastq")

&nbsp;   trimmed = trimReads(samples)

}

```



\*\*Questions to test your understanding\*\*:

1\. What does `${sample\_id}` do in the output path?

2\. What does the triple-quoted `"""..."""` section contain?

3\. How is `${reads}` different from `$reads`?



<details>

<summary>Answers</summary>



1\. \*\*String interpolation\*\*: Inserts the value of the `sample\_id` variable into the filename. If `sample\_id = "patient001"`, the output file is `"patient001\_trimmed.fastq"`



2\. \*\*Shell script\*\*: The triple-quoted section is a multi-line shell script that will be executed. It contains the actual trimmomatic command.



3\. \*\*Safety\*\*: `${reads}` clearly marks where the variable starts and ends. While `$reads` might work in simple cases, `${reads}` is safer and more explicit, especially when the variable is followed by more text or when accessing properties like `${reads.name}`



\*\*Key insight\*\*: You can now read Nextflow processes! You understand:

\- The `input:` declares what data comes in

\- The `output:` declares what's produced

\- The `script:` contains shell commands with Groovy interpolation

\- `${variable}` injects Groovy values into the script

</details>



---



\## 📝 Key Takeaways



Before moving to Day 3, make sure you understand:



✅ \*\*String interpolation with `${variable}`\*\* - You'll use this constantly  

✅ \*\*Lists with `collect {}` and `findAll {}`\*\* - Core data manipulation  

✅ \*\*The `it` keyword\*\* - Implicit parameter in closures  

✅ \*\*Maps with dot and bracket notation\*\* - For configuration and metadata  

✅ \*\*Triple quotes `"""`\*\* - For multi-line scripts in Nextflow  

✅ \*\*Groovy is similar to Python\*\* - Same logic, slightly different syntax  



\### The Mental Model



Think of Groovy as Python's cousin who:

\- Uses `${}` instead of `f"{}"`

\- Says `collect` instead of list comprehension

\- Says `findAll` instead of filter

\- Likes both `map.key` and `map\["key"]`

\- Uses `it` as a shortcut in closures

\- Is less picky about semicolons and parentheses



\### What You Actually Need to Remember



\*\*The Big 3 Patterns\*\*:



1\. \*\*String interpolation\*\*:

&nbsp;  ```groovy

&nbsp;  "${variable}"  // Always safe

&nbsp;  ```



2\. \*\*List transformation\*\*:

&nbsp;  ```groovy

&nbsp;  list.collect { it.transform() }  // Like map

&nbsp;  list.findAll { it.condition }    // Like filter

&nbsp;  ```



3\. \*\*Map access\*\*:

&nbsp;  ```groovy

&nbsp;  config.key        // Dot notation

&nbsp;  config\["key"]     // Bracket notation

&nbsp;  ```



Everything else you can look up as needed!



---



\## 🎯 Ready for Day 3?



Tomorrow, you'll write your \*\*first Nextflow process\*\*! You now have all the Groovy knowledge you need to understand process syntax.



\### Quick Preview - You'll Understand This Tomorrow!



```groovy

process sayHello {

&nbsp;   input:

&nbsp;   val name

&nbsp;   

&nbsp;   output:

&nbsp;   path "greeting.txt"

&nbsp;   

&nbsp;   script:

&nbsp;   """

&nbsp;   echo "Hello, ${name}!" > greeting.txt

&nbsp;   """

}

```



You can already identify:

\- `val name` - input parameter

\- `${name}` - string interpolation

\- `"""..."""` - multi-line script block

\- `path "greeting.txt"` - output file



\### Groovy Quick Reference Card



Save this for reference:



```groovy

// Variables

name = "Alice"

count = 42



// String interpolation

message = "Hello, ${name}"

filename = "${id}\_processed.bam"



// Lists

samples = \["s1", "s2", "s3"]

samples << "s4"                    // Add item

transformed = samples.collect { it.toUpperCase() }

filtered = samples.findAll { it.contains("1") }

samples.each { println it }



// Maps

config = \[input: "data/", threads: 8]

input = config.input               // Dot notation

input = config\["input"]            // Bracket notation



// Closures

double = { it \* 2 }

result = \[1,2,3].collect { it \* 2 }  // \[2,4,6]



// Control flow

if (x > 5) {

&nbsp;   println "Big"

} else {

&nbsp;   println "Small"

}



for (item in list) {

&nbsp;   println item

}



// You're ready for Nextflow!

```



---



\## ✅ Day 2 Completion Checklist



Before marking Day 2 complete, ensure you can:



\- \[ ] Write string interpolation with `${}`

\- \[ ] Use `collect` to transform lists

\- \[ ] Use `findAll` to filter lists

\- \[ ] Access map values with dot notation

\- \[ ] Understand the `it` keyword

\- \[ ] Read a simple Nextflow process

\- \[ ] Know when to use triple quotes

\- \[ ] Translate simple Python to Groovy



\*\*Completed Day 2?\*\* Update your `PROGRESS.md`! Tomorrow you'll write real Nextflow code! 🎉



\*\*Your progress\*\*: 2/28 days (7.1%) complete



---



\## 🔍 Bonus: Common Gotchas



Watch out for these when starting:



1\. \*\*Quotes Matter\*\*:

&nbsp;  ```groovy

&nbsp;  'No interpolation: ${var}'   // Literal string

&nbsp;  "Yes interpolation: ${var}"  // Variable replaced

&nbsp;  ```



2\. \*\*Braces Are Your Friend\*\*:

&nbsp;  ```groovy

&nbsp;  "${id}\_file.txt"   // ✅ Correct

&nbsp;  "$id\_file.txt"     // ❌ Looks for variable "id\_file"

&nbsp;  ```



3\. \*\*List Methods Are Different\*\*:

&nbsp;  ```groovy

&nbsp;  len(list)       // ❌ Python way

&nbsp;  list.size()     // ✅ Groovy way

&nbsp;  ```



4\. \*\*True/False Are Lowercase\*\*:

&nbsp;  ```groovy

&nbsp;  True   // ❌ Python

&nbsp;  true   // ✅ Groovy

&nbsp;  ```



5\. \*\*Dot vs Bracket for Maps\*\*:

&nbsp;  ```groovy

&nbsp;  map.simple\_key     // ✅ Works

&nbsp;  map."key-with-dash"  // ✅ Dot with quotes for special chars

&nbsp;  map\["any\_key"]     // ✅ Always works

&nbsp;  ```



---



\*Tomorrow: Day 3 - Your First Nextflow Process\*



\*\*You're making great progress! See you tomorrow! 🚀\*\*

