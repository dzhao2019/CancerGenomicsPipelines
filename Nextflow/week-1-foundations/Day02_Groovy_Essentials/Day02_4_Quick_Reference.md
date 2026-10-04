# Day 2: Groovy Quick Reference Card

**Print this or bookmark it for quick lookup while learning Nextflow!**

---

## 🔤 String Interpolation (MOST IMPORTANT!)

```groovy
// Double quotes = interpolation happens
name = "Alice"
message = "Hello ${name}"          // "Hello Alice"
count = 10
result = "Count: ${count * 2}"     // "Count: 20"

// Single quotes = NO interpolation
literal = 'Hello ${name}'          // Literal: "Hello ${name}"

// Use in Nextflow scripts
script:
    """
    echo "Processing ${sample}"    // Groovy interpolation
    echo "Path: \$PATH"            // Bash variable (escaped)
    """
```

**Rule:** Always use `${}` for variables, `$""` for strings.

---

## 📋 Data Structures Quick Reference

### Lists
```groovy
list = [1, 2, 3, 4, 5]
list[0]                            // 1 (first element)
list[-1]                           // 5 (last element)
list.size()                        // 5
list.isEmpty()                     // false

// Adding to list
list << 6                          // Append
list.add(7)                        // Add
list += [8, 9]                     // Add multiple

// List methods
list.collect { it * 2 }            // Transform each
list.findAll { it > 2 }            // Filter
list.sum()                         // Add all
list.first()                       // First element
list.last()                        // Last element
```

### Maps
```groovy
map = [name: "Alice", age: 30]
map.name                           // "Alice" (dot notation)
map["age"]                         // 30 (bracket notation)
map.size()                         // 2
map.isEmpty()                      // false

// Adding to map
map.city = "New York"              // Add new key
map["country"] = "USA"             // Alternative

// Map methods
map.keySet()                       // [name, age, city, country]
map.values()                       // [Alice, 30, New York, USA]
map.each { key, val -> ... }       // Iterate
map.collect { k, v -> ... }        // Transform
```

### Tuples (Lists as Key-Value Pairs)
```groovy
// In Nextflow, tuples group related data
tuple = ["sample1", "/path/sample1.fastq"]
id = tuple[0]                      // "sample1"
file = tuple[1]                    // "/path/sample1.fastq"

// Destructure in closures
files.map { (id, path) -> ... }
```

---

## 🔄 Closures and Functional Operations

```groovy
// Closure syntax
{ x -> x * 2 }                     // Lambda-like closure
{ it * 2 }                         // Using 'it' (default param)
{ a, b -> a + b }                  // Multiple parameters

// List operations
list.each { item -> ... }          // For each item
list.map { item -> ... }           // Transform each
list.collect { item -> ... }       // Same as map
list.findAll { item -> ... }       // Filter (keep matching)
list.find { item -> ... }          // Find first match
list.any { item -> ... }           // Any matches?
list.all { item -> ... }           // All match?
list.sum()                         // Sum all

// String operations
"text".toUpperCase()               // UPPERCASE
"text".toLowerCase()               // lowercase
"text".capitalize()                // Capitalize
"text".size()                      // Length
"text".split(",")                  // Split into list

// File operations
file.name                          // Filename
file.baseName                      // Name without extension
file.extension                     // File extension
file.size()                        // File size in bytes
file.exists()                      // File exists?
```

---

## 📌 Control Flow

```groovy
// If statements
if (x > 5) {
    println "Big"
} else if (x > 2) {
    println "Medium"
} else {
    println "Small"
}

// For loops
for (i in 1..5) { println i }      // 1 to 5
for (item in list) { ... }         // Over list
for (i = 0; i < 5; i++) { ... }   // C-style

// While loops
while (count < 5) { count++ }

// Switch
switch (value) {
    case "a": ...
    case "b": ...
    default: ...
}

// Try-catch
try {
    doSomething()
} catch (Exception e) {
    println "Error: ${e.message}"
}
```

---

## 🎯 Common Groovy Patterns

### Pattern 1: Elvis Operator (Default Values)
```groovy
value = params.input ?: "default"  // Use default if null
name = person?.name ?: "Unknown"   // Safe access + default
```

### Pattern 2: Safe Navigation
```groovy
size = file?.size()                // Returns null if file is null
name = obj?.property?.nested       // Chain safe access
```

### Pattern 3: Method Chaining
```groovy
result = items
    .map { it.toUpperCase() }      // Transform
    .findAll { it.size() > 2 }     // Filter
    .sort()                         // Sort
```

### Pattern 4: Range
```groovy
range = (1..10)                    // 1 to 10 inclusive
range = (1..<10)                   // 1 to 10 exclusive
range = ('a'..'z')                 // a to z
```

### Pattern 5: Collecting into Tuple
```groovy
files.map { file ->
    [file.baseName, file]          // Return list (as tuple)
}
```

---

## 🔀 Python ↔ Groovy Cheat Sheet

| Python | Groovy | Notes |
|--------|--------|-------|
| `f"{var}"` | `"${var}"` | Interpolation |
| `[1, 2]` | `[1, 2]` | Lists |
| `{"a": 1}` | `[a: 1]` | Maps |
| `for x in items:` | `items.each { x -> }` | Loops |
| `lambda x: x*2` | `{ x -> x*2 }` | Closures |
| `map(func, list)` | `list.map(closure)` | Transform |
| `filter(cond, list)` | `list.findAll(closure)` | Filter |
| `dict["key"]` | `map["key"]` | Map access |
| `dict.get("key")` | `map.key` | Map access |
| `list[0]` | `list[0]` | Indexing |
| `None` | `null` | Null value |
| `True/False` | `true/false` | Boolean |
| `len(list)` | `list.size()` | Length |
| `str.upper()` | `str.toUpperCase()` | Case change |

---

## ❌ Common Mistakes and Fixes

| Mistake | Why It's Wrong | Fix |
|---------|---|---|
| `'Hello ${name}'` | Single quotes don't interpolate | `"Hello ${name}"` |
| `{"key": "value"}` | Python dict syntax | `[key: "value"]` |
| `for x in items:` | Python syntax, missing parens/braces | `for (x in items) { }` |
| `.toUpper()` | Method doesn't exist | `.toUpperCase()` |
| `list[0]` + no error check | May get null error | Use `list?.get(0)` |
| `obj.property` on null | Crashes if null | Use `obj?.property` |
| `{ item }` | Missing `->` for parameters | `{ item -> item }` |

---

## 🧬 Nextflow-Specific Groovy

### Accessing Parameters
```groovy
// Define in params block
params {
    input = "data/"
    cores = 4
}

// Access in workflow
input_dir = params.input
threads = params.cores

// Interpolate in script
script:
    """
    tool --threads ${params.cores} ${params.input}
    """
```

### Working with Files
```groovy
// Create file object
f = file("path/to/file.txt")
f.baseName                         // file
f.extension                        // txt
f.size()                           // bytes
f.exists()                         // true/false

// In channels
Channel.fromPath("*.fastq")
    .map { file -> file.baseName }
    .view()
```

### Process Metadata Access
```groovy
process MY_PROCESS {
    input:
        path reads
    
    script:
        // Access file properties
        """
        echo "Processing ${reads.baseName}"
        """
}
```

### Simple Conditionals in Script
```groovy
process ALIGN {
    input:
        val aligner
        path reads
    
    script:
        if (aligner == "bwa") {
            """
            bwa mem reference.fa ${reads}
            """
        } else if (aligner == "bowtie2") {
            """
            bowtie2 -x reference -U ${reads}
            """
        }
}
```

---

## 📋 Method Reference: Common String Methods

```groovy
str = "hello world"

str.size()                         // 11
str.length()                       // 11
str.toUpperCase()                  // HELLO WORLD
str.toLowerCase()                  // hello world
str.capitalize()                   // Hello world
str.reverse()                      // dlrow olleh
str.split(" ")                     // [hello, world]
str.contains("world")              // true
str.startsWith("hello")            // true
str.endsWith("world")              // true
str.replaceAll("l", "L")           // heLLo worLd
str.trim()                         // Remove whitespace
str.take(5)                        // "hello"
str.drop(6)                        // "world"
```

---

## 📋 Method Reference: Common List Methods

```groovy
list = [3, 1, 4, 1, 5, 9]

list.size()                        // 6
list.isEmpty()                     // false
list.contains(4)                   // true
list.first()                       // 3
list.last()                        // 9
list.head()                        // [3, 1, 4, 1, 5]
list.tail()                        // [1, 4, 1, 5, 9]
list.reverse()                     // [9, 5, 1, 4, 1, 3]
list.sort()                        // [1, 1, 3, 4, 5, 9]
list.unique()                      // [3, 1, 4, 5, 9]
list.sum()                         // 23
list.min()                         // 1
list.max()                         // 9
list.join(",")                     // "3,1,4,1,5,9"
```

---

## 🔗 Nextflow Groovy Context

```groovy
// Workflow structure
workflow {
    // Define channels
    ch_reads = Channel.fromPath("data/*.fastq")
    
    // Transform with closures
    ch_samples = ch_reads
        .map { file -> [file.baseName, file] }
        .filter { id, file -> file.size() > 100 }
    
    // Run processes
    PROCESS1(ch_samples)
}

// Process structure
process EXAMPLE {
    input:
        tuple val(id), path(file)
    output:
        tuple val(id), path("*.result")
    script:
        """
        // Groovy: variable substitution
        echo "Processing ${id} from ${file}"
        tool ${file} > ${id}.result
        """
}
```

---

## 🎯 Most Important Rules

1. **Use `${}` for interpolation** in double-quoted strings
2. **Single quotes** = no interpolation (literals)
3. **Maps use** `[key: value]` **syntax** (not `{"key": value}`)
4. **Closures are** `{ param -> body }` or `{ it }`
5. **Methods need** `()` even with no arguments
6. **Dot notation works for maps:** `map.key` is same as `map["key"]`
7. **Elvis operator** `?: "default"` for null defaults
8. **Safe navigation** `?.` prevents null errors

---

## 🔍 Debugging Checklist

When something doesn't work:

- [ ] Am I using double quotes `"..."` for interpolation?
- [ ] Are map values using `:` not `:`? (Groovy, not Python)
- [ ] Do all methods have `()`?
- [ ] Is the closure syntax correct `{ ... }`?
- [ ] Should I use `?.` for safe access?
- [ ] Is the variable in scope?
- [ ] Did I miss a bracket or brace?

---

## 📚 Resources

- **Groovy syntax:** https://groovy-lang.org/syntax.html
- **Groovy GDK:** https://groovy-lang.org/gdk.html
- **Nextflow docs:** https://www.nextflow.io/docs/latest/

---

## ✨ Key Takeaway

**Groovy for Nextflow is simple:**

1. **String interpolation** - `${var}`
2. **Data structures** - Lists, maps, tuples
3. **Closures** - `{ param -> body }`
4. **Methods** - `.operation()`
5. **Control flow** - if/for/while

That's 80% of what you need. Everything else is details!

---

*Print this card and keep it handy while learning Nextflow!*
