# Day 2: Groovy Essentials - Index & Navigation

**Learning Date:** _______________  
**Time Spent:** _____ minutes  
**Status:** ⬜ Not Started | 🟡 In Progress | ✅ Complete

---

## 📂 Day 2 Materials

This Day 2 package contains everything you need to master Groovy for Nextflow.

### Files in This Package

1. **day-02-lesson.md** (Main Lesson)
   - Complete 30-minute learning experience
   - 10 key Groovy concepts explained
   - Real Nextflow examples throughout
   - Python comparisons for everything
   - 5 common gotchas and how to avoid them

2. **day-02-exercises.md** (Hands-On Practice)
   - 7 progressive exercises
   - From string interpolation to real processes
   - Detailed solutions with explanations
   - 5 minutes per exercise

3. **day-02-quick-reference.md** (Cheat Sheet)
   - One-page reference (print-friendly)
   - All syntax in quick lookup format
   - Common mistakes and fixes
   - Nextflow-specific Groovy patterns

4. **day-02-index.md** (This File)
   - Navigation and organization
   - Learning objectives
   - Suggested reading paths
   - Progress tracking

---

## 🎯 Learning Objectives

By the end of Day 2, you should be able to:

1. **Use string interpolation** ✓
   - Write `"Hello ${name}"` correctly
   - Understand `${}` vs single quotes
   - Interpolate expressions and method calls

2. **Work with lists and maps** ✓
   - Create and access lists
   - Create and access maps
   - Use dot notation vs bracket notation

3. **Understand and use closures** ✓
   - Read closure syntax `{ ... }`
   - Use `.collect()` to transform
   - Use `.findAll()` to filter
   - Understand the `it` keyword

4. **Read Nextflow code** ✓
   - Understand processes
   - Know how variables are substituted
   - Recognize closure patterns in channels

5. **Debug Groovy errors** ✓
   - Spot common syntax mistakes
   - Understand error messages
   - Know how to fix them

6. **Spot key differences from Python** ✓
   - String interpolation syntax
   - Map syntax
   - Loop syntax (closures vs for)
   - Method names

---

## ⏰ Suggested Reading Schedule (30 minutes)

**Time** | **Activity** | **File** | **Duration**
---|---|---|---
0-5 min | Introduction & Objectives | This file | 5 min
5-10 min | Warm-up: Python to Groovy | day-02-lesson.md | 5 min
10-20 min | Key Concepts (Parts 1-5) | day-02-lesson.md | 10 min
20-25 min | Key Concepts (Parts 6-10) | day-02-lesson.md | 5 min
25-27 min | Real Nextflow Examples | day-02-lesson.md | 2 min
27-30 min | Try first exercises | day-02-exercises.md | 3 min

**After Learning:**
- 15 min: Complete all 7 exercises
- 10 min: Check solutions
- 5 min: Review quick reference

---

## 📚 Content Breakdown

### Lesson Content (20 min of reading)

**Part 1: Introduction** (2 min)
- Why Groovy?
- Learning objectives

**Part 2: Python vs Groovy** (3 min)
- Side-by-side comparison table
- Key differences
- What you need to know

**Part 3: String Interpolation** (3 min) - MOST IMPORTANT
- Basic syntax `"${var}"`
- Expressions and methods
- Examples you'll see in Nextflow
- Critical concept: when interpolation happens

**Part 4: Variables & Basic Types** (2 min)
- Strings, numbers, booleans
- No type declarations
- Dynamic typing like Python

**Part 5: Lists** (3 min)
- Creating and accessing
- Common operations
- Methods: `.collect()`, `.findAll()`, `.each()`
- Python equivalents

**Part 6: Maps** (3 min)
- Creating and accessing
- Dot notation vs bracket notation
- Common operations
- Iteration over maps

**Part 7: Closures** (3 min)
- Closure syntax
- Using closures with lists
- The `it` parameter
- Closure examples in Nextflow

**Part 8: Groovy-specific Features** (2 min)
- Elvis operator `?:`
- Safe navigation `?.`
- Ranges `1..10`

**Part 9: Methods and Functions** (2 min)
- Defining methods
- Default parameters
- Return values

**Part 10: Control Flow** (2 min)
- If/else, for, while, switch
- Try-catch
- Differences from Python

**Real Examples** (3 min)
- 5 complete Nextflow examples
- Understanding processes
- Parameters and values
- Channels with closures

---

## 🔑 Key Concepts at a Glance

### The Three Most Important Things

1. **String Interpolation**
   ```groovy
   name = "Alice"
   message = "Hello ${name}"  // "Hello Alice"
   ```
   - Use double quotes `"..."`
   - Use `${}` for variables
   - Works with expressions too

2. **Closures**
   ```groovy
   items.collect { it * 2 }           // Transform
   items.findAll { it > 5 }           // Filter
   items.each { item -> println item }// Iterate
   ```
   - Blocks of code in `{ ... }`
   - `it` is the default parameter
   - Used everywhere in Nextflow

3. **Maps**
   ```groovy
   config = [name: "Alice", age: 30]
   config.name    // "Alice"
   config.age     // 30
   ```
   - Use `[key: value]` syntax
   - Access with `.key` or `["key"]`
   - Common in Nextflow metadata

---

## 💻 Exercise Breakdown

### Exercise 1: String Interpolation (5 min)
- Part A: Basic interpolation
- Part B: Interpolation with expressions
- Part C: Interpolation with methods
- **Skills:** Understanding `${}`

### Exercise 2: Lists and Closures (5 min)
- Part A: Using `.collect()` to transform
- Part B: Using `.findAll()` to filter
- Part C: Using `it` in closures
- **Skills:** List transformation and filtering

### Exercise 3: Maps (5 min)
- Part A: Creating and accessing maps
- Part B: Iterating over maps
- Part C: Transforming list of maps
- **Skills:** Working with map structures

### Exercise 4: Reading Nextflow Code (5 min)
- Part A: Simple process reading
- Part B: Processes with parameters
- Part C: Channel operations
- **Skills:** Understanding Nextflow patterns

### Exercise 5: Debugging Groovy Errors (5 min)
- Part A: Missing parentheses
- Part B: Single vs double quotes
- Part C: Map syntax differences
- **Skills:** Spotting and fixing errors

### Exercise 6: Writing Simple Groovy (5 min)
- Part A: Writing closures
- Part B: String interpolation
- Part C: Creating maps
- **Skills:** Writing Groovy code

### Exercise 7: Real Nextflow Pattern (10 min)
- Challenge: Understanding a complete process
- Four comprehension questions
- Real-world context
- **Skills:** Applying all concepts together

---

## ✅ Completion Checklist

**Understanding Key Concepts:**
- [ ] String interpolation: `"${variable}"`
- [ ] Lists: creating, accessing, methods
- [ ] Maps: creating, accessing with `.` and `[]`
- [ ] Closures: syntax and parameter usage
- [ ] Single vs double quotes (literals vs interpolation)
- [ ] Elvis operator: `value ?: default`
- [ ] Safe navigation: `obj?.property`

**Practice Skills:**
- [ ] Can write a closure to transform data
- [ ] Can write a closure to filter data
- [ ] Can create and access a map
- [ ] Can interpolate values into strings
- [ ] Can read simple Nextflow processes
- [ ] Can spot common Groovy mistakes

**Nextflow Context:**
- [ ] Understand how Groovy variables get substituted into scripts
- [ ] Know when interpolation happens (before script runs)
- [ ] Recognize closure patterns in channel operations
- [ ] Understand tuple syntax for grouping data
- [ ] Know how to access file properties (`.baseName`, `.size()`)

**Ready for Day 3:**
- [ ] Comfortable reading Groovy code
- [ ] Understand basic Groovy syntax
- [ ] Can explain Groovy to someone else
- [ ] Ready to write first Nextflow process

---

## 🤔 Reflection Questions

Take 5 minutes to reflect on these:

### Question 1: Comparison
How is Groovy string interpolation different from Python?

**My answer:**
```
_________________________________________________________________
_________________________________________________________________
```

### Question 2: Closures
Why might `.collect()` with a closure be useful in bioinformatics?

**My answer:**
```
_________________________________________________________________
_________________________________________________________________
```

### Question 3: Nextflow Context
When you see `${sample_id}` in a Nextflow script block, what's happening?

**My answer:**
```
_________________________________________________________________
_________________________________________________________________
```

### Question 4: Common Error
Someone writes: `config = {"name": "Alice"}` and gets an error. Why?

**My answer:**
```
_________________________________________________________________
_________________________________________________________________
```

---

## 🎯 Most Important Concepts (Priority Order)

1. ⭐⭐⭐ **String Interpolation** - Used constantly, must understand
2. ⭐⭐⭐ **Closures** - Core to all list/channel operations
3. ⭐⭐⭐ **Maps** - Used for storing metadata
4. ⭐⭐ **Lists** - Fundamental data structure
5. ⭐⭐ **Groovy Syntax** - Different from Python
6. ⭐ **Advanced Features** - Elvis operator, safe navigation (nice to have)

---

## 🔄 Comparison: Python ↔ Groovy

When you're confused, remember these key differences:

| Python | Groovy | Remember |
|--------|--------|----------|
| `f"Hello {name}"` | `"Hello ${name}"` | Use `${}` in Groovy |
| `{"key": value}` | `[key: value]` | No quotes around keys |
| `for x in items:` | `items.each { x -> }` | Groovy uses closures |
| `list[0]` | `list[0]` | Same indexing |
| `None` | `null` | Different keyword |
| `lambda x: x*2` | `{ x -> x*2 }` | Closures, not lambdas |

---

## 🚀 What's Next (Day 3 Preview)

Tomorrow you'll:
- Write your first Nextflow process
- Use everything you learned today
- Create input/output declarations
- Write script blocks with interpolation
- Understand how Groovy and bash interact

**Why Day 2 matters for Day 3:**
- Process definitions are mostly Groovy
- String interpolation creates bash commands
- Closures aren't used much in single processes, but you'll understand the syntax
- Maps organize input/output metadata

---

## 📋 Progress Tracking

**Date Started:** _______________  
**Estimated Completion:** 30-45 minutes  
**Time Actually Spent:** _____ minutes  

**Which parts did you struggle with?**
```
_________________________________________________________________
```

**Which parts clicked quickly?**
```
_________________________________________________________________
```

**What do you want to review?**
```
_________________________________________________________________
```

**Confidence Level (1-10):** _____

---

## 💬 Quick Help Guide

**"I don't understand string interpolation"**
→ Read Part 3 of lesson, then try Exercise 1

**"Closures are confusing"**
→ Read Part 7 of lesson, then try Exercise 2

**"I can't read Nextflow code"**
→ Try Exercise 4 first, then the lesson

**"I keep making syntax errors"**
→ Try Exercise 5, then check the cheat sheet

**"I need practice"**
→ Do all 7 exercises, check solutions, try variations

---

## 🎓 How to Know You're Ready for Day 3

You're ready when you can:

✅ Explain what `${var}` does  
✅ Write `.collect { it * 2 }`  
✅ Create a map and access its values  
✅ Read a simple Nextflow process  
✅ Spot differences from Python  
✅ Fix a simple Groovy syntax error  
✅ Understand closures in context  

If you can do all these, you're ready for Day 3!

---

## 🎉 Celebration Point!

You've completed Day 1 (Nextflow concepts) and now you're finishing Day 2 (Groovy syntax). **You're 2/28 days through!** 🎊

You now understand:
- What Nextflow is and why it matters (Day 1)
- How to read and write Groovy code (Day 2)
- You're prepared to write your first process (Day 3)

This is real progress!

---

## 📞 If You Get Stuck

**Can't understand a concept?**
1. Re-read the relevant section
2. Try the exercise for that concept
3. Check the solution explanation
4. Look at the cheat sheet
5. Try variations of the exercise

**Keep a notebook**
- Write down patterns you see
- Write simple practice code
- Note questions to revisit later

---

## 🔗 Resources

- **Lesson:** day-02-lesson.md
- **Exercises:** day-02-exercises.md
- **Cheat Sheet:** day-02-quick-reference.md
- **Groovy Docs:** https://groovy-lang.org/
- **Nextflow Docs:** https://www.nextflow.io/docs/latest/

---

## ✨ Key Takeaway

You don't need to be a Groovy expert. You only need to understand the small subset used in Nextflow.

**That subset:**
1. String interpolation
2. Lists and maps
3. Closures
4. Basic control flow
5. Method calls

Everything else is optional!

---

**Status: Day 2 Complete or In Progress? Check your checklist above!**

*Day 2 of 28 - You're building the foundation for writing production-ready Nextflow workflows!*
