# 🧬 Day 6: Running Workflows and Understanding Execution
**Week 1: Foundations · Saturday · Nextflow Mastery Course**

> Run pipelines confidently, understand the `work/` directory, and use `-resume` and execution reports

---

## 🗂️ File Structure

```
Day06_Running_Workflows_Execution/
├── Day06_1_Lesson.md              (Main lesson — START HERE)
├── Day06_2_InClass_Exercises.md   (In-class exercises)
├── Day06_3_Homework.md            (Homework)
├── Day06_4_Quick_Reference.md     (Cheat sheet)
├── Day06_5_Extended_Lesson.md     (Extended lesson, optional)
└── README.md                      (This file)
```

---

## ⚡ Quick Start

1. Read the **[lesson](Day06_1_Lesson.md)** (30 min)
2. Work through the **[in-class exercises](Day06_2_InClass_Exercises.md)** (~25 min)
3. Complete the **[homework](Day06_3_Homework.md)** (~45 min)
4. Keep the **[quick reference](Day06_4_Quick_Reference.md)** open while coding
5. Log the day in [`Progress_Log.md`](../../00_Course_Guide/Progress_Log.md)

---

## 📚 What Each File Contains

| File | Purpose | Time | Best For |
|------|---------|------|----------|
| [Day06_1_Lesson](Day06_1_Lesson.md) | Main lesson — **START HERE** | 30 min | Learning the concepts |
| [Day06_2_InClass_Exercises](Day06_2_InClass_Exercises.md) | Guided exercises with collapsible solutions | 25 min | Practice |
| [Day06_3_Homework](Day06_3_Homework.md) | Self-paced homework with solutions | 45 min | Consolidation |
| [Day06_4_Quick_Reference](Day06_4_Quick_Reference.md) | One-page cheat sheet | 5 min | Quick lookup while coding |
| [Day06_5_Extended_Lesson](Day06_5_Extended_Lesson.md) | Long-form version of the lesson with extra examples and built-in exercises | optional | Deeper dive / review |

---

## 🎯 Learning Objectives

By the end of Day 6 you will be able to:

- Run a pipeline with `nextflow run` and read the live console output
- Navigate a task's work directory and explain each hidden `.command.*` file
- Use `-resume` and predict which tasks will be cached
- Use `nextflow log` to find past runs and failed tasks
- Generate `-with-report`, `-with-trace`, `-with-timeline` and `-with-dag` outputs

---

## ✅ Completion Checklist

- [ ] I can read the console progress lines (run name, hash, tag, counts)
- [ ] I can find a task's work directory and read `.command.sh` / `.command.err`
- [ ] I can predict what `-resume` will re-run after a change
- [ ] I know what breaks resume (moving launch dir, deleting `work/`, touching inputs)
- [ ] I can list past runs and failed tasks with `nextflow log`
- [ ] I can generate report, trace, timeline and DAG files

---

## 🚦 Navigation

⬅️ [Day 5: Connecting Processes into Workflows](../Day05_Connecting_Processes_Workflows/README.md) · [📋 Course index](../../README.md) · [Day 7: Week 1 Review and Integration Project](../Day07_Week1_Review_Integration/README.md) ➡️

