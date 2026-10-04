# 🧬 Day 4: Understanding Channels
**Week 1: Foundations · Thursday · Nextflow Mastery Course**

> Understand channels as streams of data, create them with channel factories, and see how they drive automatic parallelism

---

## 🗂️ File Structure

```
Day04_Channels/
├── Day04_1_Lesson.md              (Main lesson — START HERE)
├── Day04_2_InClass_Exercises.md   (In-class exercises)
├── Day04_3_Homework.md            (Homework)
├── Day04_4_Quick_Reference.md     (Cheat sheet)
├── Day04_5_Extended_Lesson.md     (Extended lesson, optional)
└── README.md                      (This file)
```

---

## ⚡ Quick Start

1. Read the **[lesson](Day04_1_Lesson.md)** (30 min)
2. Work through the **[in-class exercises](Day04_2_InClass_Exercises.md)** (~25 min)
3. Complete the **[homework](Day04_3_Homework.md)** (~45 min)
4. Keep the **[quick reference](Day04_4_Quick_Reference.md)** open while coding
5. Log the day in [`Progress_Log.md`](../../00_Course_Guide/Progress_Log.md)

---

## 📚 What Each File Contains

| File | Purpose | Time | Best For |
|------|---------|------|----------|
| [Day04_1_Lesson](Day04_1_Lesson.md) | Main lesson — **START HERE** | 30 min | Learning the concepts |
| [Day04_2_InClass_Exercises](Day04_2_InClass_Exercises.md) | Guided exercises with collapsible solutions | 25 min | Practice |
| [Day04_3_Homework](Day04_3_Homework.md) | Self-paced homework with solutions | 45 min | Consolidation |
| [Day04_4_Quick_Reference](Day04_4_Quick_Reference.md) | One-page cheat sheet | 5 min | Quick lookup while coding |
| [Day04_5_Extended_Lesson](Day04_5_Extended_Lesson.md) | Long-form version of the lesson with extra examples and built-in exercises | optional | Deeper dive / review |

---

## 🎯 Learning Objectives

By the end of Day 4 you will be able to:

- Explain a channel as a **stream** (conveyor belt), not a list
- Distinguish **queue channels** (consumed once) from **value channels** (reused forever)
- Create channels with `Channel.of`, `Channel.fromPath`, `Channel.fromFilePairs`, `Channel.value`
- Inspect channel contents with `.view()`
- Predict how many tasks a process will run given its input channels

---

## ✅ Completion Checklist

- [ ] I can explain why a channel is a stream, not a list
- [ ] I know the difference between queue and value channels
- [ ] I can explain (and fix) the "only one task ran" bug
- [ ] I can create channels with `of`, `fromPath`, `fromFilePairs`, `value`
- [ ] I use `.view()` to inspect channel contents
- [ ] I can predict the number of tasks from channel sizes

---

## 🚦 Navigation

⬅️ [Day 3: Your First Nextflow Process](../Day03_First_Process/README.md) · [📋 Course index](../../README.md) · [Day 5: Connecting Processes into Workflows](../Day05_Connecting_Processes_Workflows/README.md) ➡️

