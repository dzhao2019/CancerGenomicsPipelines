# Day 3: Your First Nextflow Process - Index & Navigation

**Learning Date:** _______________  
**Time Spent:** _____ minutes  
**Status:** ⬜ Not Started | 🟡 In Progress | ✅ Complete

---

## 📂 Day 3 Materials

Everything you need to write your first Nextflow process.

### Files in This Package

1. **day-03-lesson.md** (Main Lesson)
   - Complete 30-minute learning experience
   - 10 key concepts about processes
   - Real bioinformatics examples
   - Detailed process anatomy
   - Common error patterns

2. **day-03-exercises.md** (Hands-On Practice)
   - 6 progressive exercises
   - From understanding structure to debugging
   - Detailed solutions with explanations
   - Real Nextflow patterns

3. **day-03-quick-reference.md** (Cheat Sheet)
   - Print-friendly one-page reference
   - Process template (copy & modify)
   - Common patterns
   - Error fixes
   - Input/output types

4. **day-03-index.md** (This File)
   - Navigation and organization
   - Learning objectives
   - Suggested reading schedule
   - Progress tracking

---

## 🎯 Learning Objectives

By the end of Day 3, you should be able to:

1. **Understand process anatomy** ✓
   - 5 parts: name, directives, input, output, script
   - What each part does
   - Why each part exists

2. **Write input declarations** ✓
   - Know when to use `path` vs `val`
   - Understand tuples for grouped data
   - Use correct syntax

3. **Write output declarations** ✓
   - Use glob patterns correctly
   - Match patterns to created files
   - Preserve metadata through tuples

4. **Write script blocks** ✓
   - Use Groovy variable interpolation
   - Chain commands together
   - Use special variables like `${task.cpus}`

5. **Run a process** ✓
   - Understand the work directory
   - Know what happens at each step
   - See your code in action

6. **Debug process errors** ✓
   - Identify syntax errors
   - Fix output mismatches
   - Resolve input problems

7. **Write production-ready code** ✓
   - Use best practices
   - Name things clearly
   - Organize outputs

---

## ⏰ Suggested Reading Schedule (30 minutes)

**Time** | **Activity** | **File** | **Duration**
---|---|---|---
0-5 min | Introduction & Objectives | This file | 5 min
5-12 min | Process Structure (Concepts 1-5) | day-03-lesson.md | 7 min
12-20 min | Examples and Patterns (Concepts 6-10) | day-03-lesson.md | 8 min
20-27 min | Try Exercises 1-3 | day-03-exercises.md | 7 min
27-30 min | Quick Reference Overview | day-03-quick-reference.md | 3 min

**After Learning:**
- 15 min: Complete Exercises 4-6
- 10 min: Check solutions
- 5 min: Review quick reference

---

## 📚 Content Breakdown

### Lesson: 10 Key Concepts

1. **Process Structure: The Anatomy** (5 parts)
   - Process declaration
   - Directives (optional)
   - Input block
   - Output block
   - Script block

2. **Understanding Inputs: Data Types**
   - `path` - files
   - `val` - simple values
   - `tuple` - grouped data
   - `stdin` - pipes

3. **Understanding Outputs: What You Create**
   - `path` with patterns
   - `val` for values
   - `stdout` for output
   - Pattern matching rules

4. **The Script Block: Where Work Happens**
   - Variable interpolation
   - Multi-line commands
   - Language mixing
   - Bash/Python/R

5. **Example: Simple Quality Control Process**
   - Real FastQC example
   - Step-by-step breakdown
   - Key directives

6. **Process Communication: How Data Flows**
   - Input channels
   - Output channels
   - Data passing

7. **Real Example: Alignment Process**
   - Complete real example
   - Multiple tools (BWA, samtools)
   - Best practices

8. **Common Process Patterns**
   - QC pattern
   - Transformation pattern
   - Analysis pattern
   - Aggregation pattern

9. **Understanding the Work Directory**
   - Where processes run
   - Why it matters
   - Debugging with work dir

10. **Debugging: Common Errors**
    - File not found
    - Syntax errors
    - Input missing
    - Tool not installed

---

## ✅ Completion Checklist

**Understanding Key Concepts:**
- [ ] 5-part process structure
- [ ] When to use each directive
- [ ] `path` vs `val` inputs
- [ ] Glob patterns for outputs
- [ ] Groovy variable interpolation
- [ ] Work directory concept
- [ ] Tuple grouping
- [ ] Dynamic file naming

**Practice Skills:**
- [ ] Can fix syntax errors
- [ ] Can fix output mismatches
- [ ] Can write simple process
- [ ] Can name outputs dynamically
- [ ] Can preserve sample IDs
- [ ] Can use special variables
- [ ] Can debug process errors

**Ready for Day 4:**
- [ ] Comfortable with process structure
- [ ] Understand all 5 parts
- [ ] Can write a working process
- [ ] Know common patterns
- [ ] Ready to connect processes

---

## 🔑 Key Concepts at a Glance

### The 5 Parts of a Process
1. **Directives** - Tell Nextflow resources needed
2. **Input** - Declare what you need (`path` for files, `val` for values)
3. **Output** - Declare what you create (patterns like `*.bam`)
4. **Script** - Do the actual work (bash, Python, R, etc.)
5. **Name** - Describe what it does (UPPERCASE)

### Most Important Rules
- Input `path` = file, `val` = simple value
- Output pattern must match created files
- Use `${}` for Groovy variables
- Process names are UPPERCASE
- Return sample IDs in output tuples

---

## 🎯 Most Important Concepts (Priority Order)

1. ⭐⭐⭐ **Input/Output Types**
   - Know when to use `path` vs `val`
   - Output patterns must match files

2. ⭐⭐⭐ **Script Block**
   - Where actual work happens
   - Groovy interpolation

3. ⭐⭐ **Directives**
   - Tell Nextflow what's needed
   - Optional but recommended

4. ⭐⭐ **Preserving Data**
   - Use tuples to keep sample IDs
   - Matters for downstream processes

5. ⭐ **Advanced Features**
   - Special variables (task.cpus, etc.)
   - publishDir for organization

---

## 🧠 Reflection Questions

Take 5 minutes to think about these:

### Question 1: Process Purpose
Think of a tool you use (BWA, FastQC, etc.):
- What are its inputs?
- What are its outputs?
- Would you use `path` or `val` for each?

**My answer:**
```
_________________________________________________________________
_________________________________________________________________
```

### Question 2: Output Naming
Why is it important to include sample ID in output filenames?

**My answer:**
```
_________________________________________________________________
_________________________________________________________________
```

### Question 3: Error Debugging
What would you do if a process fails with "No files match pattern"?

**My answer:**
```
_________________________________________________________________
_________________________________________________________________
```

---

## 📖 Suggested Reading Order

1. **START:** This file (day-03-index.md) - 5 min
2. **MAIN:** day-03-lesson.md - 20 min
3. **PRACTICE:** Exercises 1-3 - 7 min
4. **PRACTICE:** Exercises 4-6 - 15 min
5. **REFERENCE:** day-03-quick-reference.md - 3 min

---

## 🚀 What's Next (Day 4)

Tomorrow you'll learn about **Channels** - connecting processes together.

**You already have:**
- Understanding of process structure ✅
- Ability to write processes ✅
- Knowledge of inputs/outputs ✅

**Tomorrow you'll learn:**
- How to create channels
- How to pass data between processes
- How to build workflows

**Your first complete workflow is tomorrow!**

---

## 💻 Code You'll Write Today

By end of Day 3, you'll write:
- At least 3 working processes
- Handle inputs/outputs correctly
- Use dynamic naming
- Preserve sample IDs
- Format for production use

---

## 🎓 How You Know You're Ready for Day 4

✅ Can write `process PROCESS_NAME {}`  
✅ Understand `path` and `val` inputs  
✅ Can write output patterns  
✅ Can interpolate variables in scripts  
✅ Know how directives work  
✅ Can debug common errors  
✅ Feel confident about processes  

---

## 📊 Progress Tracking

**Date Started:** _______________  
**Time Spent:** _____ minutes  
**Exercises Completed:** ___/6  

**What felt easy:**
```
_________________________________________________________________
```

**What felt challenging:**
```
_________________________________________________________________
```

**I want to review:**
```
_________________________________________________________________
```

**Confidence Level (1-10):** _____

---

## 🎉 Milestone

You're on **Day 3 of 28** (11% through course!)

**Week 1 Progress:**
- Day 1: What Nextflow is ✅
- Day 2: Groovy syntax ✅
- Day 3: Your first process (TODAY)
- Days 4-5: Workflows and channels
- Days 6-7: Practical applications

**By end of Week 1:** You'll have built complete workflows!

---

## 💡 Remember

A process is:
- **Not** a function in Python (though similar)
- **Not** magic - just a declared interface
- **Isolated** - can run on any computer
- **Reusable** - write once, use for many inputs
- **Parallelizable** - Nextflow runs them efficiently

---

## 🔗 Resources

- **Lesson:** day-03-lesson.md
- **Exercises:** day-03-exercises.md
- **Quick ref:** day-03-quick-reference.md
- **Official docs:** https://www.nextflow.io/docs/latest/process.html

---

## ✨ Key Takeaway

**A Nextflow process is:**
1. A declaration of what you need (inputs)
2. A specification of what you create (outputs)
3. Some code that does the work (script)
4. Metadata about resources (directives)
5. All wrapped in Nextflow's orchestration magic

Master this pattern, and you can write any bioinformatics workflow!

---

**Status: Ready to begin Day 3? Check your understanding of Nextflow from Days 1-2, then dive into the lesson!**

*Day 3 of 28 - You're about to write your first Nextflow code! 🚀*

