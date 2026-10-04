# 📋 Nextflow Mastery: Course Outline

> A 28-day Nextflow DSL2 curriculum for Python programmers — from "why not just a Python loop?" to a tested, scalable, monitored RNA-seq pipeline.

This outline matches the course folder exactly. Every day title links to that day's `README.md`; the file links open the materials directly. For the compact day list see the [master index](../README.md).

---

## Course Overview

| | |
|---|---|
| **Duration** | 28 days (4 weeks, 7 days per week) |
| **Audience** | Bioinformaticians and Python programmers new to workflow managers |
| **Prerequisites** | Python basics; familiarity with FASTQ/BAM/VCF and the command line |
| **Core path** | 30-minute lesson per day → 14 hours |
| **Full path** | Lesson 30 min + in-class exercises ~25 min + homework ~45 min ≈ 1.5–2 h per day → ~45 hours |
| **Outcome** | Design, build, test, deploy and operate production Nextflow pipelines |

### Every day has the same materials

| File | Purpose | Time |
|---|---|---|
| `README.md` | Day overview, links, objectives, checklist | 2 min |
| `DayNN_1_Lesson.md` | Main lesson — start here | 30 min |
| `DayNN_2_InClass_Exercises.md` | Guided exercises (recognition → modification → creation), solutions collapsed | ~25 min |
| `DayNN_3_Homework.md` | Self-paced homework with solutions and reflection questions | ~45 min |
| `DayNN_4_Quick_Reference.md` | One-page cheat sheet | — |
| `DayNN_5_Extended_Lesson.md` | Long-form lesson with extra examples (Days 1–14, 16, 17) | optional |
| `DayNN_extra_*.md` | Supplementary material (Days 1, 17, 21) | optional |

### Teaching approach

- **Python bridge** — every concept is mapped to a Python equivalent, including where the analogy breaks
- **Real bioinformatics** — FastQC, Trim Galore, BWA, SAMtools, GATK, STAR, featureCounts, DESeq2, MultiQC, nf-core modules
- **Progressive projects** — one pipeline grows from Week 1 to Week 3; Week 4 builds a new one to production standard
- **Central idea** — *the data-dependency graph is the program*

---

## Course at a Glance

| Day | Weekday | Topic | Milestone |
|---|---|---|---|
| 1 | Mon | [What Nextflow Actually Is](../Week1_Foundations/Day01_What_Nextflow_Is/README.md) |  |
| 2 | Tue | [Groovy Essentials for Nextflow](../Week1_Foundations/Day02_Groovy_Essentials/README.md) |  |
| 3 | Wed | [Your First Nextflow Process](../Week1_Foundations/Day03_First_Process/README.md) |  |
| 4 | Thu | [Understanding Channels](../Week1_Foundations/Day04_Channels/README.md) |  |
| 5 | Fri | [Connecting Processes into Workflows](../Week1_Foundations/Day05_Connecting_Processes_Workflows/README.md) |  |
| 6 | Sat | [Running Workflows and Understanding Execution](../Week1_Foundations/Day06_Running_Workflows_Execution/README.md) |  |
| 7 | Sun | [Week 1 Review and Integration Project](../Week1_Foundations/Day07_Week1_Review_Integration/README.md) | 🏗️ Week 1 project |
| 8 | Mon | [Handling Parameters and Making Workflows Flexible](../Week2_Practical_Pipelines/Day08_Parameters_Flexibility/README.md) |  |
| 9 | Tue | [Working with Collections and Channel Operators](../Week2_Practical_Pipelines/Day09_Channel_Operators/README.md) |  |
| 10 | Wed | [Combining Multiple Inputs in Processes](../Week2_Practical_Pipelines/Day10_Combining_Multiple_Inputs/README.md) |  |
| 11 | Thu | [Publishing and Managing Outputs](../Week2_Practical_Pipelines/Day11_Publishing_Outputs/README.md) |  |
| 12 | Fri | [Error Handling and Robustness](../Week2_Practical_Pipelines/Day12_Error_Handling/README.md) |  |
| 13 | Sat | [Working with Containers for Reproducibility](../Week2_Practical_Pipelines/Day13_Containers/README.md) |  |
| 14 | Sun | [Debugging Workflows Systematically](../Week2_Practical_Pipelines/Day14_Debugging/README.md) | 🏗️ Week 2 project |
| 15 | Mon | [Subworkflows and Modularity](../Week3_Advanced_Patterns/Day15_Subworkflows_Modularity/README.md) |  |
| 16 | Tue | [Importing and Using Existing Workflows](../Week3_Advanced_Patterns/Day16_Importing_nf-core/README.md) |  |
| 17 | Wed | [Advanced Channel Operations](../Week3_Advanced_Patterns/Day17_Advanced_Channel_Operations/README.md) |  |
| 18 | Thu | [Conditional Execution and Control Flow](../Week3_Advanced_Patterns/Day18_Conditional_Execution/README.md) |  |
| 19 | Fri | [Performance Optimization and Resource Management](../Week3_Advanced_Patterns/Day19_Performance_Resources/README.md) |  |
| 20 | Sat | [Configuration Files and Profiles](../Week3_Advanced_Patterns/Day20_Configuration_Profiles/README.md) |  |
| 21 | Sun | [Sharing and Publishing Workflows](../Week3_Advanced_Patterns/Day21_Sharing_Publishing/README.md) | 🏗️ Week 3 project (`minivar`) |
| 22 | Mon | [Building a Complete RNA-seq Pipeline — Part 1: Design, Input and QC](../Week4_Production/Day22_RNAseq_Part1_QC/README.md) | 🏗️ `rnaseq-mini` Part 1 |
| 23 | Tue | [Building the RNA-seq Pipeline — Part 2: Alignment and Quantification](../Week4_Production/Day23_RNAseq_Part2_Align_Quant/README.md) | 🏗️ `rnaseq-mini` Part 2 |
| 24 | Wed | [Building the RNA-seq Pipeline — Part 3: Differential Expression and Reporting](../Week4_Production/Day24_RNAseq_Part3_DiffExp/README.md) | 🏗️ `rnaseq-mini` Part 3 — end to end |
| 25 | Thu | [Testing, Validation and Quality Assurance](../Week4_Production/Day25_Testing_Validation/README.md) | 🏗️ `rnaseq-mini` tested |
| 26 | Fri | [Scaling to Cluster and Cloud Environments](../Week4_Production/Day26_Scaling_Cluster_Cloud/README.md) | 🏗️ `rnaseq-mini` scaled |
| 27 | Sat | [Monitoring, Logging and Troubleshooting in Production](../Week4_Production/Day27_Monitoring_Production/README.md) | 🏗️ `rnaseq-mini` monitored |
| 28 | Sun | [Final Review and Continuing Your Journey](../Week4_Production/Day28_Final_Review/README.md) | 🎓 Portfolio |

---

## 🗓️ Week 1: Foundations

**Folder**: [`Week1_Foundations/`](../Week1_Foundations/)  
**Theme**: {theme}

> **Week 1 project (Day 7):** a multi-sample QC pipeline — FastQC, read counting with PASS/FAIL, trimming of passing samples, summary table and MultiQC.

### [Day 1: What Nextflow Actually Is](../Week1_Foundations/Day01_What_Nextflow_Is/README.md)
**Monday** · [Lesson](../Week1_Foundations/Day01_What_Nextflow_Is/Day01_1_Lesson.md) · [In-class](../Week1_Foundations/Day01_What_Nextflow_Is/Day01_2_InClass_Exercises.md) · [Homework](../Week1_Foundations/Day01_What_Nextflow_Is/Day01_3_Homework.md) · [Quick ref](../Week1_Foundations/Day01_What_Nextflow_Is/Day01_4_Quick_Reference.md) · [Extended](../Week1_Foundations/Day01_What_Nextflow_Is/Day01_5_Extended_Lesson.md) · [Extended Lesson Solutions](../Week1_Foundations/Day01_What_Nextflow_Is/Day01_extra_Extended_Lesson_Solutions.md) · [Groovy Examples](../Week1_Foundations/Day01_What_Nextflow_Is/Day01_extra_Groovy_Examples.md)

**Goal**: Understand what Nextflow is, which problem it solves, and when to choose it over a Python script

**Learning objectives**:
- Distinguish **data processing** (what a tool does) from **workflow orchestration** (how tools are coordinated)
- Name the three building blocks of every Nextflow pipeline: **processes, channels, workflows**
- Explain why Nextflow parallelises automatically while a Python `for` loop does not
- Describe **resumability**, **portability** and **reproducibility** in one sentence each
- Decide whether a given task is better solved with Python, Nextflow, or both

**Lesson covers**: The Pipeline Coordination Problem · The Three Building Blocks · The Four Superpowers

**In-class**: Python, Nextflow, or Both? · Identify the Building Blocks · Compute the Speed-up · Sketch Your Own Pipeline

**Homework**: Translate a Python Script into Building Blocks · Declarative Thinking · Make the Case

---

### [Day 2: Groovy Essentials for Nextflow](../Week1_Foundations/Day02_Groovy_Essentials/README.md)
**Tuesday** · [Lesson](../Week1_Foundations/Day02_Groovy_Essentials/Day02_1_Lesson.md) · [In-class](../Week1_Foundations/Day02_Groovy_Essentials/Day02_2_InClass_Exercises.md) · [Homework](../Week1_Foundations/Day02_Groovy_Essentials/Day02_3_Homework.md) · [Quick ref](../Week1_Foundations/Day02_Groovy_Essentials/Day02_4_Quick_Reference.md) · [Extended](../Week1_Foundations/Day02_Groovy_Essentials/Day02_5_Extended_Lesson.md)

**Goal**: Learn just enough Groovy to read and write Nextflow scripts — mapped directly onto Python you already know

**Learning objectives**:
- Use **string interpolation** (`"${sample}.bam"`) and know when single quotes disable it
- Create and manipulate **lists** and **maps** (Python lists and dicts)
- Write **closures** (`{ x -> x * 2 }`, `{ it * 2 }`) — Groovy's lambdas
- Use `collect`, `findAll`, `each` the way you use comprehensions and `filter` in Python
- Recognise Groovy's optional parentheses, which make Nextflow code look "magic"

**Lesson covers**: Why Groovy? · Variables and Strings · Lists, Maps and Closures

**In-class**: String Interpolation · Lists and Closures · Maps · Reading Nextflow Code · Debugging Groovy Errors · Writing Simple Groovy · Real Nextflow Pattern

**Homework**: Sample Sheet Wrangling · Fix the Quoting Bugs · Parse File Names into Tuples

---

### [Day 3: Your First Nextflow Process](../Week1_Foundations/Day03_First_Process/README.md)
**Wednesday** · [Lesson](../Week1_Foundations/Day03_First_Process/Day03_1_Lesson.md) · [In-class](../Week1_Foundations/Day03_First_Process/Day03_2_InClass_Exercises.md) · [Homework](../Week1_Foundations/Day03_First_Process/Day03_3_Homework.md) · [Quick ref](../Week1_Foundations/Day03_First_Process/Day03_4_Quick_Reference.md) · [Extended](../Week1_Foundations/Day03_First_Process/Day03_5_Extended_Lesson.md)

**Goal**: Write, run and understand a complete Nextflow process that wraps a real bioinformatics tool

**Learning objectives**:
- Name the parts of a process: **directives**, `input:`, `output:`, `script:`
- Choose the right input qualifier: `val`, `path`, `tuple`
- Declare outputs so Nextflow can find the files a tool produced
- Run a single-process workflow and find its results in `work/`
- Explain why a process is *not* a Python function

**Lesson covers**: Process Anatomy · Inputs and Outputs · The Script Block and Running a Process

**In-class**: Understanding Process Structure · Fix the Broken Process · Write a Simple Process · Reading Real Nextflow Code · Debugging Real Errors · Your First Real Process

**Homework**: Wrap FastQC · Fix Five Broken Processes · Python Inside a Process

---

### [Day 4: Understanding Channels](../Week1_Foundations/Day04_Channels/README.md)
**Thursday** · [Lesson](../Week1_Foundations/Day04_Channels/Day04_1_Lesson.md) · [In-class](../Week1_Foundations/Day04_Channels/Day04_2_InClass_Exercises.md) · [Homework](../Week1_Foundations/Day04_Channels/Day04_3_Homework.md) · [Quick ref](../Week1_Foundations/Day04_Channels/Day04_4_Quick_Reference.md) · [Extended](../Week1_Foundations/Day04_Channels/Day04_5_Extended_Lesson.md)

**Goal**: Understand channels as streams of data, create them with channel factories, and see how they drive automatic parallelism

**Learning objectives**:
- Explain a channel as a **stream** (conveyor belt), not a list
- Distinguish **queue channels** (consumed once) from **value channels** (reused forever)
- Create channels with `Channel.of`, `Channel.fromPath`, `Channel.fromFilePairs`, `Channel.value`
- Inspect channel contents with `.view()`
- Predict how many tasks a process will run given its input channels

**Lesson covers**: What Is a Channel? · Queue vs Value Channels · Channel Factories

**In-class**: Predict the Output · Count the Tasks · Build the Right Channels · Debug: "Why did only one sample run?"

**Homework**: Channel Explorer Script · Queue vs Value Experiment · Choose the Factory

---

### [Day 5: Connecting Processes into Workflows](../Week1_Foundations/Day05_Connecting_Processes_Workflows/README.md)
**Friday** · [Lesson](../Week1_Foundations/Day05_Connecting_Processes_Workflows/Day05_1_Lesson.md) · [In-class](../Week1_Foundations/Day05_Connecting_Processes_Workflows/Day05_2_InClass_Exercises.md) · [Homework](../Week1_Foundations/Day05_Connecting_Processes_Workflows/Day05_3_Homework.md) · [Quick ref](../Week1_Foundations/Day05_Connecting_Processes_Workflows/Day05_4_Quick_Reference.md) · [Extended](../Week1_Foundations/Day05_Connecting_Processes_Workflows/Day05_5_Extended_Lesson.md)

**Goal**: Chain processes into a multi-step pipeline by passing outputs as inputs

**Learning objectives**:
- Use a process's **output channel** as the input of the next process
- Access outputs with assignment, `.out`, and named `emit:` outputs
- Use the pipe operator `|` for simple linear chains
- Fan out (one channel → several processes) and fan in (`.collect()`, `.mix()`)
- Read a workflow block as a **dependency graph**, not a sequence of calls

**Lesson covers**: Outputs Are Channels · Accessing Outputs · Fan-Out and Fan-In

**In-class**: Draw the Graph · Fix the Shape Mismatches · Add a Step · Design: Week 1 Mini Pipeline

**Homework**: Three-Process QC Chain · Named Outputs Refactor · Fan-Out, Fan-In Design

---

### [Day 6: Running Workflows and Understanding Execution](../Week1_Foundations/Day06_Running_Workflows_Execution/README.md)
**Saturday** · [Lesson](../Week1_Foundations/Day06_Running_Workflows_Execution/Day06_1_Lesson.md) · [In-class](../Week1_Foundations/Day06_Running_Workflows_Execution/Day06_2_InClass_Exercises.md) · [Homework](../Week1_Foundations/Day06_Running_Workflows_Execution/Day06_3_Homework.md) · [Quick ref](../Week1_Foundations/Day06_Running_Workflows_Execution/Day06_4_Quick_Reference.md) · [Extended](../Week1_Foundations/Day06_Running_Workflows_Execution/Day06_5_Extended_Lesson.md)

**Goal**: Run pipelines confidently, understand the `work/` directory, and use `-resume` and execution reports

**Learning objectives**:
- Run a pipeline with `nextflow run` and read the live console output
- Navigate a task's work directory and explain each hidden `.command.*` file
- Use `-resume` and predict which tasks will be cached
- Use `nextflow log` to find past runs and failed tasks
- Generate `-with-report`, `-with-trace`, `-with-timeline` and `-with-dag` outputs

**Lesson covers**: Running a Pipeline · The Work Directory · Resume · Logs and Reports (bonus reading)

**In-class**: Read the Console · Predict the Resume · Investigate a Failure · Generate and Interpret Reports

**Homework**: Break It, Fix It, Resume It · Work Directory Scavenger Hunt · Execution Reports

---

### [Day 7: Week 1 Review and Integration Project](../Week1_Foundations/Day07_Week1_Review_Integration/README.md)
**Sunday** · [Lesson](../Week1_Foundations/Day07_Week1_Review_Integration/Day07_1_Lesson.md) · [In-class](../Week1_Foundations/Day07_Week1_Review_Integration/Day07_2_InClass_Exercises.md) · [Homework](../Week1_Foundations/Day07_Week1_Review_Integration/Day07_3_Homework.md) · [Quick ref](../Week1_Foundations/Day07_Week1_Review_Integration/Day07_4_Quick_Reference.md) · [Extended](../Week1_Foundations/Day07_Week1_Review_Integration/Day07_5_Extended_Lesson.md)

**Goal**: Consolidate Week 1 by reviewing every core concept and assembling them into one working multi-sample QC pipeline

**Milestone**: 🏗️ Week 1 project

**Learning objectives**:
- Explain the Week 1 mental model in one diagram: **processes + channels + workflow**
- Recall the key syntax from each day without looking it up
- Design a multi-sample QC pipeline from a written specification
- Run, break, and resume that pipeline, and inspect its work directories
- Identify your own weak spots before Week 2

**Lesson covers**: Week 1 in One Page · Integration Project Specification · Reference Implementation Walkthrough

**In-class**: Build the Skeleton · Add Status, Trimming and Summary · Break, Inspect, Resume · Code Review

**Homework**: Finish and Harden the QC Pipeline · Explain It Back · Plan the Week 2 Upgrade

---


## 🔧 Week 2: Practical Pipelines

**Folder**: [`Week2_Practical_Pipelines/`](../Week2_Practical_Pipelines/)  
**Theme**: {theme}

> **Week 2 project (Day 14):** a robust QC → align → call-variants pipeline with parameters, containers, error handling, organised outputs and a diagnostic report.

### [Day 8: Handling Parameters and Making Workflows Flexible](../Week2_Practical_Pipelines/Day08_Parameters_Flexibility/README.md)
**Monday** · [Lesson](../Week2_Practical_Pipelines/Day08_Parameters_Flexibility/Day08_1_Lesson.md) · [In-class](../Week2_Practical_Pipelines/Day08_Parameters_Flexibility/Day08_2_InClass_Exercises.md) · [Homework](../Week2_Practical_Pipelines/Day08_Parameters_Flexibility/Day08_3_Homework.md) · [Quick ref](../Week2_Practical_Pipelines/Day08_Parameters_Flexibility/Day08_4_Quick_Reference.md) · [Extended](../Week2_Practical_Pipelines/Day08_Parameters_Flexibility/Day08_5_Extended_Lesson.md)

**Goal**: Transform hard-coded workflows into configurable, reusable pipelines

**Learning objectives**:
- Declare `params` with sensible defaults
- Pass parameters on the command line without editing code
- Access `params` inside processes and channel factories
- Validate required parameters and fail fast with clear error messages
- Explain the Python analogy: `params` ≈ `argparse` + config dict

**Lesson covers**: Why Parameters? · Core `params` Syntax · Validation Patterns

**In-class**: Spot the Problems · Add Parameters to an Existing Workflow · Design the Parameter Block

**Homework**: Build a Configurable QC Pipeline · Parameter Inheritance via Config File · Robust Validation Challenge

---

### [Day 9: Working with Collections and Channel Operators](../Week2_Practical_Pipelines/Day09_Channel_Operators/README.md)
**Tuesday** · [Lesson](../Week2_Practical_Pipelines/Day09_Channel_Operators/Day09_1_Lesson.md) · [In-class](../Week2_Practical_Pipelines/Day09_Channel_Operators/Day09_2_InClass_Exercises.md) · [Homework](../Week2_Practical_Pipelines/Day09_Channel_Operators/Day09_3_Homework.md) · [Quick ref](../Week2_Practical_Pipelines/Day09_Channel_Operators/Day09_4_Quick_Reference.md) · [Extended](../Week2_Practical_Pipelines/Day09_Channel_Operators/Day09_5_Extended_Lesson.md)

**Goal**: Master the operator toolkit that turns raw channels into precisely shaped data streams

**Learning objectives**:
- Use `map {}` to transform every element in a channel
- Use `filter {}` to selectively pass elements downstream
- Use `groupTuple` to collect related elements by key
- Use `flatten` and `toList` / `collect` to reshape channel structure
- Chain multiple operators into readable pipelines
- Explain why operators are lazily evaluated (and why that matters)

**Lesson covers**: The Operator Mental Model · Core Operators · Operator Chaining · Complete Worked Example

**In-class**: Read the Chain, Predict the Output · Fix the Broken Chain · Build a Sample-Grouping Pipeline

**Homework**: Samplesheet-Driven Pipeline · Multi-Condition Grouping · Operator Chain Debugging

---

### [Day 10: Combining Multiple Inputs in Processes](../Week2_Practical_Pipelines/Day10_Combining_Multiple_Inputs/README.md)
**Wednesday** · [Lesson](../Week2_Practical_Pipelines/Day10_Combining_Multiple_Inputs/Day10_1_Lesson.md) · [In-class](../Week2_Practical_Pipelines/Day10_Combining_Multiple_Inputs/Day10_2_InClass_Exercises.md) · [Homework](../Week2_Practical_Pipelines/Day10_Combining_Multiple_Inputs/Day10_3_Homework.md) · [Quick ref](../Week2_Practical_Pipelines/Day10_Combining_Multiple_Inputs/Day10_4_Quick_Reference.md) · [Extended](../Week2_Practical_Pipelines/Day10_Combining_Multiple_Inputs/Day10_5_Extended_Lesson.md)

**Goal**: Master every pattern for getting multiple pieces of data into a single process — the foundation of every real alignment and variant-calling pipeline

**Learning objectives**:
- Use `tuple` inputs to keep sample metadata paired with file paths across process boundaries
- Pass a shared reference file to every sample using a value channel
- Use `each` to run a process over every combination of samples × conditions
- Wire multi-input processes correctly in the `workflow {}` block
- Carry tuples through a multi-step pipeline without losing the sample ID

**Lesson covers**: Tuple inputs: keep sample ID + file together · Shared references via value channels · `each` for combinatorial inputs · Worked example: BWA → SAMtools → GATK

**In-class**: Trace the Tuple · Fix the Broken Pipeline · Design and Wire a Two-Tool Comparison

**Homework**: The Full Alignment Pipeline · The `each` Benchmarking Pipeline · Multi-Output Process Routing

---

### [Day 11: Publishing and Managing Outputs](../Week2_Practical_Pipelines/Day11_Publishing_Outputs/README.md)
**Thursday** · [Lesson](../Week2_Practical_Pipelines/Day11_Publishing_Outputs/Day11_1_Lesson.md) · [In-class](../Week2_Practical_Pipelines/Day11_Publishing_Outputs/Day11_2_InClass_Exercises.md) · [Homework](../Week2_Practical_Pipelines/Day11_Publishing_Outputs/Day11_3_Homework.md) · [Quick ref](../Week2_Practical_Pipelines/Day11_Publishing_Outputs/Day11_4_Quick_Reference.md) · [Extended](../Week2_Practical_Pipelines/Day11_Publishing_Outputs/Day11_5_Extended_Lesson.md)

**Goal**: Transform a working pipeline into a professional one by mastering how Nextflow saves, organises, and exposes results

**Learning objectives**:
- Add `publishDir` to any process to save outputs to a user-facing directory
- Choose the right publish mode: `copy`, `symlink`, `move`, or `link`
- Use `pattern:` to publish only a subset of a process's outputs
- Apply multiple `publishDir` directives to route different file types to different folders
- Explain the difference between intermediate outputs (work directory only) and final results (published)
- Design a logical, navigable output directory structure for a multi-step pipeline

**Lesson covers**: The Work Directory Problem · `publishDir` Fundamentals · Selective Publishing with `pattern:` · Multiple `publishDir` Blocks · Output Directory Design

**In-class**: Publishing Audit · Fix and Extend Publishing · Design an Output Structure

**Homework**: Retrofit Publishing onto a Complete Pipeline · Per-Sample Directory Structure · The `saveAs` Renaming Challenge

---

### [Day 12: Error Handling and Robustness](../Week2_Practical_Pipelines/Day12_Error_Handling/README.md)
**Friday** · [Lesson](../Week2_Practical_Pipelines/Day12_Error_Handling/Day12_1_Lesson.md) · [In-class](../Week2_Practical_Pipelines/Day12_Error_Handling/Day12_2_InClass_Exercises.md) · [Homework](../Week2_Practical_Pipelines/Day12_Error_Handling/Day12_3_Homework.md) · [Quick ref](../Week2_Practical_Pipelines/Day12_Error_Handling/Day12_4_Quick_Reference.md) · [Extended](../Week2_Practical_Pipelines/Day12_Error_Handling/Day12_5_Extended_Lesson.md)

**Goal**: Write pipelines that survive real-world failures — transient cluster errors, bad samples, resource exhaustion — without crashing the entire run

**Learning objectives**:
- Apply `errorStrategy` to control what happens when a process task fails
- Configure `retry` with `maxRetries` to handle transient failures automatically
- Write dynamic `errorStrategy` closures that respond to exit codes
- Use `ignore` and `finish` to skip failed samples without stopping the pipeline
- Escalate resources (memory, CPUs) on retry to defeat OOM failures
- Distinguish between permanent errors (bad data) and transient errors (cluster hiccup)

**Lesson covers**: Why Pipelines Fail in Practice · The Four `errorStrategy` Values · Dynamic Strategies with Closures · Resource Escalation on Retry · Global Error Strategy in Config

**In-class**: Failure Taxonomy · Fix the Fragile Pipeline · Design a Full Error Policy

**Homework**: Retrofit Error Handling onto a Complete Pipeline · Global Config Error Policy · Failure Logging and Reporting

---

### [Day 13: Working with Containers for Reproducibility](../Week2_Practical_Pipelines/Day13_Containers/README.md)
**Saturday** · [Lesson](../Week2_Practical_Pipelines/Day13_Containers/Day13_1_Lesson.md) · [In-class](../Week2_Practical_Pipelines/Day13_Containers/Day13_2_InClass_Exercises.md) · [Homework](../Week2_Practical_Pipelines/Day13_Containers/Day13_3_Homework.md) · [Quick ref](../Week2_Practical_Pipelines/Day13_Containers/Day13_4_Quick_Reference.md) · [Extended](../Week2_Practical_Pipelines/Day13_Containers/Day13_5_Extended_Lesson.md)

**Goal**: Make every process in your pipeline reproducible on any system — laptop, HPC, or cloud — by specifying exactly which software it needs

**Learning objectives**:
- Explain why containers solve the "works on my machine" problem in bioinformatics
- Add a `container` directive to any process
- Find the right image for a bioinformatics tool using BioContainers and Quay.io
- Enable Docker or Singularity in `nextflow.config`
- Use different containers for different processes in the same pipeline
- Apply the nf-core container convention for consistent image naming

**Lesson covers**: The Reproducibility Problem · Container Concepts · The `container` Directive · Finding the Right Image · Enabling Containers in Config

**In-class**: Spot the Container Problems · Fix the Container Configuration · Containerise a Complete Pipeline

**Homework**: Fully Containerise the Week 2 Pipeline · RNA-seq Pipeline with Per-Process Containers · Container Lookup Research Exercise

---

### [Day 14: Debugging Workflows Systematically](../Week2_Practical_Pipelines/Day14_Debugging/README.md)
**Sunday** · [Lesson](../Week2_Practical_Pipelines/Day14_Debugging/Day14_1_Lesson.md) · [In-class](../Week2_Practical_Pipelines/Day14_Debugging/Day14_2_InClass_Exercises.md) · [Homework](../Week2_Practical_Pipelines/Day14_Debugging/Day14_3_Homework.md) · [Quick ref](../Week2_Practical_Pipelines/Day14_Debugging/Day14_4_Quick_Reference.md) · [Extended](../Week2_Practical_Pipelines/Day14_Debugging/Day14_5_Extended_Lesson.md)

**Goal**: Develop a repeatable, fast debugging protocol that turns "the pipeline failed" into a diagnosed and fixed problem within minutes

**Milestone**: 🏗️ Week 2 project

**Learning objectives**:
- Read `.nextflow.log` and locate the exact failure point
- Navigate to a failed task's work directory and inspect its hidden files
- Reconstruct and manually reproduce any failed command
- Use `-with-trace` to find performance bottlenecks
- Apply the simplification strategy — test with small data, stub downstream processes
- Diagnose the 8 most common Nextflow failure patterns by symptom

**Lesson covers**: The Debugging Mindset · The 5-Step Investigation Protocol · Anatomy of the Work Directory · Reading `.nextflow.log` · Simplification Strategies · The 8 Common Failure Patterns

**In-class**: Crime Scene Investigation · Read and Diagnose a Trace File · Fix a Systematically Broken Pipeline

**Homework**: Full Diagnostic Report · Build a Debuggable Pipeline · The Week 2 Capstone Audit

---


## 🚀 Week 3: Advanced Patterns

**Folder**: [`Week3_Advanced_Patterns/`](../Week3_Advanced_Patterns/)  
**Theme**: {theme}

> **Week 3 project (Day 21):** `minivar` — the Week 2 pipeline refactored into a modular, schema-validated, versioned GitHub repository anyone can run with `-r v1.0.0`.

### [Day 15: Subworkflows and Modularity](../Week3_Advanced_Patterns/Day15_Subworkflows_Modularity/README.md)
**Monday** · [Lesson](../Week3_Advanced_Patterns/Day15_Subworkflows_Modularity/Day15_1_Lesson.md) · [In-class](../Week3_Advanced_Patterns/Day15_Subworkflows_Modularity/Day15_2_InClass_Exercises.md) · [Homework](../Week3_Advanced_Patterns/Day15_Subworkflows_Modularity/Day15_3_Homework.md) · [Quick ref](../Week3_Advanced_Patterns/Day15_Subworkflows_Modularity/Day15_4_Quick_Reference.md)

**Goal**: Break a monolithic `main.nf` into reusable, testable, composable pieces

**Learning objectives**:
- Define a named workflow with `take:`, `main:`, and `emit:` blocks
- Import processes and subworkflows across files with `include { } from`
- Alias an imported component so it can be invoked more than once
- Lay out a project using `modules/local/` and `subworkflows/local/`
- Test a subworkflow in isolation with `nextflow run main.nf -entry NAME`

**Lesson covers**: The anonymous workflow vs. named workflows · The three blocks · Calling a subworkflow and reading its outputs · `include` — the import statement · Aliasing — the rule that surprises everyone · Project layout · Testing a subworkflow alone

**In-class**: Recognition: Trace the Modular Project · Modification: Extract a Subworkflow · Creation: Build an Aliased Subworkflow

**Homework**: Refactor Your Week 2 Pipeline · Nested Subworkflows and Aliased Composition · Debug the Broken Modular Project

---

### [Day 16: Importing and Using Existing Workflows](../Week3_Advanced_Patterns/Day16_Importing_nf-core/README.md)
**Tuesday** · [Lesson](../Week3_Advanced_Patterns/Day16_Importing_nf-core/Day16_1_Lesson.md) · [In-class](../Week3_Advanced_Patterns/Day16_Importing_nf-core/Day16_2_InClass_Exercises.md) · [Homework](../Week3_Advanced_Patterns/Day16_Importing_nf-core/Day16_3_Homework.md) · [Quick ref](../Week3_Advanced_Patterns/Day16_Importing_nf-core/Day16_4_Quick_Reference.md) · [Extended](../Week3_Advanced_Patterns/Day16_Importing_nf-core/Day16_5_Extended_Lesson.md)

**Goal**: Stop writing modules by hand and start using the ~1,500 the community already maintains

**Learning objectives**:
- Run an existing nf-core pipeline with version pinning and test profiles
- Install, update, and remove nf-core modules with `nf-core modules`
- Read an nf-core module and explain every block in it
- Adapt your channels to the `tuple val(meta), path(files)` convention
- Configure a module with `ext.args` **without editing the module file**
- Explain what `modules.json` locks and why that matters

**Lesson covers**: Two ways to use nf-core · The nf-core tools CLI · Anatomy of an nf-core module · The meta map — the pattern that makes it all compose · `ext.args` — configure without editing · `modules.json` — the lockfile · Versions: an ecosystem mid-migration

**In-class**: Recognition: Dissect a Real nf-core Module · Modification: Adapt Your Pipeline to nf-core Modules · Creation: Configure Without Editing

**Homework**: Replace Your Local Modules with nf-core Modules · Read a Production Pipeline · Debug the Broken nf-core Integration

---

### [Day 17: Advanced Channel Operations](../Week3_Advanced_Patterns/Day17_Advanced_Channel_Operations/README.md)
**Wednesday** · [Lesson](../Week3_Advanced_Patterns/Day17_Advanced_Channel_Operations/Day17_1_Lesson.md) · [In-class](../Week3_Advanced_Patterns/Day17_Advanced_Channel_Operations/Day17_2_InClass_Exercises.md) · [Homework](../Week3_Advanced_Patterns/Day17_Advanced_Channel_Operations/Day17_3_Homework.md) · [Quick ref](../Week3_Advanced_Patterns/Day17_Advanced_Channel_Operations/Day17_4_Quick_Reference.md) · [Extended](../Week3_Advanced_Patterns/Day17_Advanced_Channel_Operations/Day17_5_Extended_Lesson.md) · [Channel Patterns Exercises](../Week3_Advanced_Patterns/Day17_Advanced_Channel_Operations/Day17_extra_Channel_Patterns_Exercises.md) · [Channel Patterns Lesson](../Week3_Advanced_Patterns/Day17_Advanced_Channel_Operations/Day17_extra_Channel_Patterns_Lesson.md) · [Channel Patterns Quick Reference](../Week3_Advanced_Patterns/Day17_Advanced_Channel_Operations/Day17_extra_Channel_Patterns_Quick_Reference.md)

**Goal**: Drive a pipeline from a samplesheet and route, join and recombine meta-map channels like nf-core pipelines do

**Learning objectives**:
- Parse a CSV samplesheet with `splitCsv(header: true)` into `[meta, files]` tuples
- Recombine channels by key with `join` (and know what happens to unmatched keys)
- Route samples into named sub-channels with `branch`, and fork one stream with `multiMap`
- Pair tumour/normal samples with `combine(by:)` / `join` on a shared key
- Use `groupKey` so `groupTuple` emits as soon as each group is complete

**Lesson covers**: Samplesheet → Meta Map · Recombining by Key · Routing

**In-class**: Predict the Join · Samplesheet to Meta Map · Route and Rejoin · Tumour/Normal Pairing

**Homework**: Samplesheet-Driven QC Subworkflow · Multi-Lane Merge with `groupKey` · Route by Data Type

---

### [Day 18: Conditional Execution and Control Flow](../Week3_Advanced_Patterns/Day18_Conditional_Execution/README.md)
**Thursday** · [Lesson](../Week3_Advanced_Patterns/Day18_Conditional_Execution/Day18_1_Lesson.md) · [In-class](../Week3_Advanced_Patterns/Day18_Conditional_Execution/Day18_2_InClass_Exercises.md) · [Homework](../Week3_Advanced_Patterns/Day18_Conditional_Execution/Day18_3_Homework.md) · [Quick ref](../Week3_Advanced_Patterns/Day18_Conditional_Execution/Day18_4_Quick_Reference.md)

**Goal**: Decide *which steps run* based on parameters and data — without breaking resume, channels or readability

**Learning objectives**:
- Skip or enable steps with `params.skip_*` flags and workflow-level `if` blocks
- Keep downstream wiring valid when a step is skipped (`Channel.empty()`, `mix`)
- Declare **optional outputs** and **optional inputs**
- Explain the `when:` directive and why nf-core prefers workflow-level `if` / `ext.when`
- Choose between parameter-driven and data-driven (`branch`) conditionals

**Lesson covers**: Parameter-Driven Steps · Optional Inputs and Outputs · `when:` and Data-Driven Conditionals

**In-class**: Graph-Time or Run-Time? · Add Skip Flags Without Breaking Wiring · Optional BED Input · Fail Fast When Everything Fails QC

**Homework**: Aligner Switch · Single-End / Paired-End Routing · Conditional Steps with Config Only

---

### [Day 19: Performance Optimization and Resource Management](../Week3_Advanced_Patterns/Day19_Performance_Resources/README.md)
**Friday** · [Lesson](../Week3_Advanced_Patterns/Day19_Performance_Resources/Day19_1_Lesson.md) · [In-class](../Week3_Advanced_Patterns/Day19_Performance_Resources/Day19_2_InClass_Exercises.md) · [Homework](../Week3_Advanced_Patterns/Day19_Performance_Resources/Day19_3_Homework.md) · [Quick ref](../Week3_Advanced_Patterns/Day19_Performance_Resources/Day19_4_Quick_Reference.md)

**Goal**: Request the right CPUs, memory and time for each process, recover automatically from out-of-memory failures, and find bottlenecks from execution reports

**Learning objectives**:
- Set `cpus`, `memory` and `time` directives and use `task.cpus` / `task.memory` in scripts
- Group resource requests with **labels** (`process_low`, `process_high`, …)
- Retry OOM-killed tasks with dynamically increasing resources (`task.attempt`)
- Cap requests to what the machine/cluster offers with `resourceLimits`
- Read the execution report and trace to find and fix bottlenecks

**Lesson covers**: Resource Directives · Labels · Dynamic Resources and Retry · Profiling and Common Optimisations

**In-class**: Spot the Resource Problems · Read the Report · Write the Retry Config · Too Many Tiny Tasks

**Homework**: Label Your Pipeline · Profile and Right-Size · Input-Size-Aware Memory

---

### [Day 20: Configuration Files and Profiles](../Week3_Advanced_Patterns/Day20_Configuration_Profiles/README.md)
**Saturday** · [Lesson](../Week3_Advanced_Patterns/Day20_Configuration_Profiles/Day20_1_Lesson.md) · [In-class](../Week3_Advanced_Patterns/Day20_Configuration_Profiles/Day20_2_InClass_Exercises.md) · [Homework](../Week3_Advanced_Patterns/Day20_Configuration_Profiles/Day20_3_Homework.md) · [Quick ref](../Week3_Advanced_Patterns/Day20_Configuration_Profiles/Day20_4_Quick_Reference.md)

**Goal**: Separate *what* the pipeline does (code) from *where and how* it runs (configuration), using layered config files and profiles

**Learning objectives**:
- Structure `nextflow.config` with scopes (`params`, `process`, `docker`, `singularity`, `executor`, …)
- Split configuration with `includeConfig` into `base.config`, `modules.config`, `test.config`
- Define and combine **profiles** (`-profile test,docker`, `-profile slurm,singularity`)
- Move `publishDir` and `ext.args` into `modules.config` with `withName` selectors — fixing the Day 15 aliasing problem
- Predict which value wins using the configuration precedence rules

**Lesson covers**: Config Anatomy · Profiles · `modules.config`: Configure Processes by Name · Precedence: Who Wins?

**In-class**: Who Wins? · Fix the Aliasing Publish Problem · Design a Profile Set · `ext.args` Instead of Code Edits

**Homework**: Refactor Into Layered Config · A Site Config for Your HPC · Params File for Reproducible Runs

---

### [Day 21: Sharing and Publishing Workflows](../Week3_Advanced_Patterns/Day21_Sharing_Publishing/README.md)
**Sunday** · [Lesson](../Week3_Advanced_Patterns/Day21_Sharing_Publishing/Day21_1_Lesson.md) · [In-class](../Week3_Advanced_Patterns/Day21_Sharing_Publishing/Day21_2_InClass_Exercises.md) · [Homework](../Week3_Advanced_Patterns/Day21_Sharing_Publishing/Day21_3_Homework.md) · [Quick ref](../Week3_Advanced_Patterns/Day21_Sharing_Publishing/Day21_4_Quick_Reference.md) · [Variant Calling Capstone Exercises](../Week3_Advanced_Patterns/Day21_Sharing_Publishing/Day21_extra_Variant_Calling_Capstone_Exercises.md) · [Variant Calling Capstone Lesson](../Week3_Advanced_Patterns/Day21_Sharing_Publishing/Day21_extra_Variant_Calling_Capstone_Lesson.md) · [Variant Calling Capstone Quick Reference](../Week3_Advanced_Patterns/Day21_Sharing_Publishing/Day21_extra_Variant_Calling_Capstone_Quick_Reference.md)

**Goal**: Package your pipeline so anyone can find, understand, run and cite a specific version of it

**Milestone**: 🏗️ Week 3 project (`minivar`)

**Learning objectives**:
- Lay out a pipeline repository the way the community expects (nf-core style)
- Write a README and `docs/` that answer "how do I run it?" and "what do I get?"
- Describe and validate parameters with `nextflow_schema.json` and the **nf-schema** plugin
- Version with Git tags, `manifest.version`, `CHANGELOG.md` and semantic versioning
- Run a shared pipeline directly from GitHub with `nextflow run owner/repo -r v1.0.0`

**Lesson covers**: Repository Layout · Documentation · Parameter Schema and Validation · Versioning and Release

**In-class**: Repo Review · Write Schema Entries · Version Bump Decisions · Run a Shared Pipeline

**Homework**: Week 3 Project: Make It Shareable · Add Schema Validation · Release v1.1.0

---


## 🏆 Week 4: Production Workflows

**Folder**: [`Week4_Production/`](../Week4_Production/)  
**Theme**: {theme}

> **Week 4 project (Days 22–27):** `rnaseq-mini` — FastQC → Trim Galore → STAR → SAMtools → featureCounts → DESeq2 → report + MultiQC, with nf-test/CI, SLURM and cloud profiles, and production monitoring.

### [Day 22: Building a Complete RNA-seq Pipeline — Part 1: Design, Input and QC](../Week4_Production/Day22_RNAseq_Part1_QC/README.md)
**Monday** · [Lesson](../Week4_Production/Day22_RNAseq_Part1_QC/Day22_1_Lesson.md) · [In-class](../Week4_Production/Day22_RNAseq_Part1_QC/Day22_2_InClass_Exercises.md) · [Homework](../Week4_Production/Day22_RNAseq_Part1_QC/Day22_3_Homework.md) · [Quick ref](../Week4_Production/Day22_RNAseq_Part1_QC/Day22_4_Quick_Reference.md)

**Goal**: Design a production RNA-seq pipeline end to end, then build its input-handling and QC/preprocessing subworkflow

**Milestone**: 🏗️ `rnaseq-mini` Part 1

**Learning objectives**:
- Design a multi-stage pipeline as **subworkflows with explicit contracts** (`take:` / `emit:` shapes)
- Set up the `rnaseq-mini` project skeleton you'll extend on Days 23–28
- Parse an RNA-seq samplesheet with `strandedness` and `condition` into meta maps
- Build a `QC_TRIM` subworkflow: FastQC (raw) → Trim Galore → FastQC (trimmed)
- Collect QC outputs for MultiQC with the accumulator pattern

**Lesson covers**: Pipeline Design · Samplesheet and Meta Map · The QC_TRIM Subworkflow

**In-class**: Write the Contracts · Samplesheet Schema · Build QC_TRIM · Predict the Task Count

**Homework**: Complete the Part 1 Pipeline · Add a Read-Count Gate · Document Stage 1

---

### [Day 23: Building the RNA-seq Pipeline — Part 2: Alignment and Quantification](../Week4_Production/Day23_RNAseq_Part2_Align_Quant/README.md)
**Tuesday** · [Lesson](../Week4_Production/Day23_RNAseq_Part2_Align_Quant/Day23_1_Lesson.md) · [In-class](../Week4_Production/Day23_RNAseq_Part2_Align_Quant/Day23_2_InClass_Exercises.md) · [Homework](../Week4_Production/Day23_RNAseq_Part2_Align_Quant/Day23_3_Homework.md) · [Quick ref](../Week4_Production/Day23_RNAseq_Part2_Align_Quant/Day23_4_Quick_Reference.md)

**Goal**: Add the `ALIGN_QUANT` subworkflow — STAR alignment, BAM indexing and statistics, and per-sample gene counts with featureCounts

**Milestone**: 🏗️ `rnaseq-mini` Part 2

**Learning objectives**:
- Handle reference inputs (FASTA, GTF, STAR index) as value channels with their own meta (`[meta2, path]`)
- Build the STAR index only when one isn't supplied — and cache it with `storeDir`
- Wire nf-core `star/align`, `samtools/index`, `samtools/stats` and `subread/featurecounts`
- Use `join` to keep BAM + BAI together and attach the GTF per sample
- Drive per-sample strandedness through `meta.strandedness`

**Lesson covers**: Reference Handling · Alignment and BAM Processing · Quantification with featureCounts

**In-class**: Find the Channel Bugs · Strandedness Detective · Build ALIGN_QUANT · Reference Parameters

**Homework**: Finish and Run Part 2 · Mixed Single/Paired Cohort · Automatic Strandedness Sanity Check

---

### [Day 24: Building the RNA-seq Pipeline — Part 3: Differential Expression and Reporting](../Week4_Production/Day24_RNAseq_Part3_DiffExp/README.md)
**Wednesday** · [Lesson](../Week4_Production/Day24_RNAseq_Part3_DiffExp/Day24_1_Lesson.md) · [In-class](../Week4_Production/Day24_RNAseq_Part3_DiffExp/Day24_2_InClass_Exercises.md) · [Homework](../Week4_Production/Day24_RNAseq_Part3_DiffExp/Day24_3_Homework.md) · [Quick ref](../Week4_Production/Day24_RNAseq_Part3_DiffExp/Day24_4_Quick_Reference.md)

**Goal**: Turn per-sample counts into differential expression results, plots and a final report — and connect the full pipeline end to end

**Milestone**: 🏗️ `rnaseq-mini` Part 3 — end to end

**Learning objectives**:
- Gather per-sample outputs into one matrix with `collect` and a Python helper in `bin/`
- Run an R/DESeq2 analysis as a process with its own container
- Pass a **contrast** and sample metadata into an analysis step cleanly
- Produce plots (PCA, MA, volcano, heatmap) and a self-contained HTML report
- Wire QC → ALIGN_QUANT → DIFFEXP → MULTIQC into the finished `RNASEQ_MINI` workflow

**Lesson covers**: Merge Counts · DESeq2 as a Process · Report · The Finished Workflow

**In-class**: Why Did MERGE_COUNTS Re-run? · Build the DIFFEXP Subworkflow · Multiple Contrasts in Parallel · Container Check

**Homework**: Run the Complete Pipeline · Report Generator · Write `docs/output.md` for Parts 2–3

---

### [Day 25: Testing, Validation and Quality Assurance](../Week4_Production/Day25_Testing_Validation/README.md)
**Thursday** · [Lesson](../Week4_Production/Day25_Testing_Validation/Day25_1_Lesson.md) · [In-class](../Week4_Production/Day25_Testing_Validation/Day25_2_InClass_Exercises.md) · [Homework](../Week4_Production/Day25_Testing_Validation/Day25_3_Homework.md) · [Quick ref](../Week4_Production/Day25_Testing_Validation/Day25_4_Quick_Reference.md)

**Goal**: Prove the pipeline works — and keeps working — with nf-test, snapshots, stub runs and continuous integration

**Milestone**: 🏗️ `rnaseq-mini` tested

**Learning objectives**:
- Explain the testing pyramid for pipelines: stub → module → subworkflow → pipeline
- Write **nf-test** tests for a process, a subworkflow and the whole pipeline
- Use **snapshots** to detect unexpected output changes, and handle non-deterministic files
- Use `stub:` blocks and `-stub-run` for second-scale wiring tests
- Run tests automatically on every push with GitHub Actions

**Lesson covers**: What to Test · nf-test Basics · Subworkflow, Pipeline and Stub Tests · Continuous Integration

**In-class**: Which Test Catches It? · Write a Process Test · Stub the Pipeline · Snapshot Review

**Homework**: Subworkflow Tests for ALIGN_QUANT · A Truth-Based DESeq2 Test · CI + Badge

---

### [Day 26: Scaling to Cluster and Cloud Environments](../Week4_Production/Day26_Scaling_Cluster_Cloud/README.md)
**Friday** · [Lesson](../Week4_Production/Day26_Scaling_Cluster_Cloud/Day26_1_Lesson.md) · [In-class](../Week4_Production/Day26_Scaling_Cluster_Cloud/Day26_2_InClass_Exercises.md) · [Homework](../Week4_Production/Day26_Scaling_Cluster_Cloud/Day26_3_Homework.md) · [Quick ref](../Week4_Production/Day26_Scaling_Cluster_Cloud/Day26_4_Quick_Reference.md)

**Goal**: Run the same tested pipeline on an HPC scheduler and on cloud batch services by changing only configuration

**Milestone**: 🏗️ `rnaseq-mini` scaled

**Learning objectives**:
- Explain how an **executor** turns tasks into jobs (local, SLURM, AWS Batch, Google Batch)
- Write a SLURM + Apptainer/Singularity profile with sensible queue and submission limits
- Write AWS Batch and Google Cloud Batch profiles with object-storage work directories
- Decide where the Nextflow **head job** should run, and how data moves (staging, Fusion, Wave)
- Estimate and control cost and throttling at scale

**Lesson covers**: Executors and the Head Job · HPC with SLURM · Cloud: AWS Batch and Google Cloud Batch

**In-class**: Diagnose the Failed HPC Launch · Write the SLURM Profile · Translate to AWS · Where Should It Run?

**Homework**: Run `rnaseq-mini` on Your HPC · Cloud Readiness Review · Cost and Throughput Estimate

---

### [Day 27: Monitoring, Logging and Troubleshooting in Production](../Week4_Production/Day27_Monitoring_Production/README.md)
**Saturday** · [Lesson](../Week4_Production/Day27_Monitoring_Production/Day27_1_Lesson.md) · [In-class](../Week4_Production/Day27_Monitoring_Production/Day27_2_InClass_Exercises.md) · [Homework](../Week4_Production/Day27_Monitoring_Production/Day27_3_Homework.md) · [Quick ref](../Week4_Production/Day27_Monitoring_Production/Day27_4_Quick_Reference.md)

**Goal**: See what a long production run is doing, get told when it finishes or fails, and diagnose failures quickly and systematically

**Milestone**: 🏗️ `rnaseq-mini` monitored

**Learning objectives**:
- Enable reports, trace, timeline and DAG permanently in config, published to `pipeline_info/`
- Write `workflow.onComplete` / `onError` handlers that summarise a run and send notifications
- Use `nextflow log` queries to find failed, slow or retried tasks across runs
- Monitor runs live with Seqera Platform
- Apply a production troubleshooting playbook to the most common failure classes

**Lesson covers**: Always-On Observability · Completion Hooks and Notifications · Forensics and Live Monitoring · Production Troubleshooting Playbook

**In-class**: Read the Trace · Write the Completion Hook · `nextflow log` Queries · Triage Five Failures

**Homework**: Production-Ready Observability for `rnaseq-mini` · Failure Drills · Write the Runbook

---

### [Day 28: Final Review and Continuing Your Journey](../Week4_Production/Day28_Final_Review/README.md)
**Sunday** · [Lesson](../Week4_Production/Day28_Final_Review/Day28_1_Lesson.md) · [In-class](../Week4_Production/Day28_Final_Review/Day28_2_InClass_Exercises.md) · [Homework](../Week4_Production/Day28_Final_Review/Day28_3_Homework.md) · [Quick ref](../Week4_Production/Day28_Final_Review/Day28_4_Quick_Reference.md)

**Goal**: Consolidate four weeks into one mental model, assess your competencies honestly, assemble your portfolio, and plan your next steps

**Milestone**: 🎓 Portfolio

**Learning objectives**:
- Explain Nextflow end to end with one central idea: *the data-dependency graph is the program*
- Map every course topic to the layer of a pipeline it belongs to
- Rate yourself against the course competencies and find your gaps
- Assemble a portfolio of the pipelines you built
- Choose concrete next steps: community, projects, deeper topics

**Lesson covers**: The Whole Course on One Page · Best Practices That Matter Most · Self-Assessment · Portfolio and Next Steps

**In-class**: Concept Map · Explain It to a Python Colleague · Code Review: Spot the Anti-Patterns · Next-Steps Pitch

**Homework**: Build Your Portfolio Document · Your 90-Day Learning Plan · Course Retrospective

---

## 🧭 Competency Checklist

Rate yourself 1–5 at the end of each week (Day 28 repeats this as the final self-assessment).

**Week 1 — Core**
- [ ] Explain Nextflow vs Python scripts and when to use each (D1)
- [ ] Read/write Groovy strings, lists, maps, closures (D2)
- [ ] Write processes with tuple I/O, escaping and stubs (D3)
- [ ] Choose queue vs value channels; predict task counts (D4)
- [ ] Chain processes with `emit`, `collect`, `mix` (D5)
- [ ] Inspect `work/`, use `-resume`, read reports (D6)

**Week 2 — Practical**
- [ ] Parameterise and validate (D8)
- [ ] Transform channels with operators (D9)
- [ ] Combine samples, references and `each` (D10)
- [ ] Publish organised outputs (D11)
- [ ] Handle errors and retries (D12)
- [ ] Containerise every process (D13)
- [ ] Debug systematically (D14)

**Week 3 — Advanced**
- [ ] Build subworkflows, `include`, aliasing (D15)
- [ ] Use nf-core modules, meta maps, `ext.args` (D16)
- [ ] Samplesheets, `join`, `branch`, `groupKey` (D17)
- [ ] Conditional steps, optional inputs/outputs (D18)
- [ ] Labels, dynamic resources, profiling (D19)
- [ ] Layered config and profiles (D20)
- [ ] Schema validation, versioning, sharing (D21)

**Week 4 — Production**
- [ ] Design a multi-stage pipeline with subworkflow contracts (D22–24)
- [ ] Test with nf-test, snapshots, stubs and CI (D25)
- [ ] Run on SLURM and cloud batch executors (D26)
- [ ] Monitor and troubleshoot production runs (D27)
- [ ] Assemble a portfolio and next-steps plan (D28)

---

## 🐍 Python → Nextflow Bridges

| Python | Nextflow | Introduced |
|---|---|---|
| function | process | Day 3 |
| list / generator | queue channel | Day 4 |
| constant / global | value channel | Day 4 |
| `main()` calling functions in order | `workflow {}` declaring a graph | Day 5 |
| checkpoint files | `-resume` | Day 6 |
| `argparse` | `params` + schema validation | Days 8, 21 |
| list comprehension / `filter` | `map` / `filter` operators | Day 9 |
| try/except + retry | `errorStrategy`, `maxRetries` | Day 12 |
| virtualenv / conda env | per-process containers | Day 13 |
| `import` | `include` (+ aliasing) | Day 15 |
| `pd.read_csv` / `merge` / `groupby` | `splitCsv` / `join` / `groupTuple` | Day 17 |
| config.yaml + env switch | `nextflow.config` + `-profile` | Day 20 |
| pytest | nf-test | Day 25 |
| `sbatch` / boto3 scripts | executors | Day 26 |
| `logging`, `atexit` | reports, trace, `onComplete` | Day 27 |

---

## 📝 Daily Progress Template

Record each day in [`Progress_Log.md`](Progress_Log.md):

```markdown
### Day [N]: [Topic]
**Date**: [ ]   **Time spent**: [ ] min

- [ ] Lesson   - [ ] In-class exercises   - [ ] Homework   - [ ] Quick reference reviewed

**Key learnings**:
- 

**Python → Nextflow insight**:
- 

**Challenges and how I solved them**:
- 

**Questions for later**:
- 
```

---

## 🏁 After the Course

Day 28 covers next steps in detail: community channels (nf-core Slack, community.seqera.io), practice projects (germline WGS vs nf-core/sarek, tumour/normal calling, contributing an nf-core module), and deeper topics (strict syntax, workflow outputs, Wave/Fusion, Nextflow plugins). Check the official docs for current features — the ecosystem moves fast.

---

*Course version: 2.0 — rewritten October 2026 to match the reorganised 28-day folder. Supersedes `OUTLINE.md` / `Course_Outline.md` v1.0 (November 2025), archived in `_archive/superseded/`.*  
*Maintained by: dzhao*
