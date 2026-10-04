# 🧬 Day 5: Connecting Processes into Workflows
**Week 1: Foundations · Friday · Nextflow Mastery Course**

> Chain processes into a multi-step pipeline by passing outputs as inputs

---

## 🗂️ File Structure

```
Day05_Connecting_Processes_Workflows/
├── Day05_1_Lesson.md              (Main lesson — START HERE)
├── Day05_2_InClass_Exercises.md   (In-class exercises)
├── Day05_3_Homework.md            (Homework)
├── Day05_4_Quick_Reference.md     (Cheat sheet)
├── Day05_5_Extended_Lesson.md     (Extended lesson, optional)
└── README.md                      (This file)
```

---

## ⚡ Quick Start

1. Read the **[lesson](Day05_1_Lesson.md)** (30 min)
2. Work through the **[in-class exercises](Day05_2_InClass_Exercises.md)** (~25 min)
3. Complete the **[homework](Day05_3_Homework.md)** (~45 min)
4. Keep the **[quick reference](Day05_4_Quick_Reference.md)** open while coding
5. Log the day in [`Progress_Log.md`](../../00_Course_Guide/Progress_Log.md)

---

## 📚 What Each File Contains

| File | Purpose | Time | Best For |
|------|---------|------|----------|
| [Day05_1_Lesson](Day05_1_Lesson.md) | Main lesson — **START HERE** | 30 min | Learning the concepts |
| [Day05_2_InClass_Exercises](Day05_2_InClass_Exercises.md) | Guided exercises with collapsible solutions | 25 min | Practice |
| [Day05_3_Homework](Day05_3_Homework.md) | Self-paced homework with solutions | 45 min | Consolidation |
| [Day05_4_Quick_Reference](Day05_4_Quick_Reference.md) | One-page cheat sheet | 5 min | Quick lookup while coding |
| [Day05_5_Extended_Lesson](Day05_5_Extended_Lesson.md) | Long-form version of the lesson with extra examples and built-in exercises | optional | Deeper dive / review |

---

## 🎯 Learning Objectives

By the end of Day 5 you will be able to:

- Use a process's **output channel** as the input of the next process
- Access outputs with assignment, `.out`, and named `emit:` outputs
- Use the pipe operator `|` for simple linear chains
- Fan out (one channel → several processes) and fan in (`.collect()`, `.mix()`)
- Read a workflow block as a **dependency graph**, not a sequence of calls

---

## ✅ Completion Checklist

- [ ] I can chain two processes by passing an output channel as an input
- [ ] I check that output tuple shape matches the next input
- [ ] I can use `.out`, `.out.name` and `emit:`
- [ ] I know when the pipe `|` works and when it doesn't
- [ ] I can fan out to several processes and fan in with `.collect()` / `.mix()`
- [ ] I can draw the dependency graph of a workflow block

---

## 🚦 Navigation

⬅️ [Day 4: Understanding Channels](../Day04_Channels/README.md) · [📋 Course index](../../README.md) · [Day 6: Running Workflows and Understanding Execution](../Day06_Running_Workflows_Execution/README.md) ➡️

