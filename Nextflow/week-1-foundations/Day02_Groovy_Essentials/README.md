# 🧬 Day 2: Groovy Essentials for Nextflow
**Week 1: Foundations · Tuesday · Nextflow Mastery Course**

> Learn just enough Groovy to read and write Nextflow scripts — mapped directly onto Python you already know

---

## 🗂️ File Structure

```
Day02_Groovy_Essentials/
├── Day02_1_Lesson.md              (Main lesson — START HERE)
├── Day02_2_InClass_Exercises.md   (In-class exercises)
├── Day02_3_Homework.md            (Homework)
├── Day02_4_Quick_Reference.md     (Cheat sheet)
├── Day02_5_Extended_Lesson.md     (Extended lesson, optional)
└── README.md                      (This file)
```

---

## ⚡ Quick Start

1. Read the **[lesson](Day02_1_Lesson.md)** (30 min)
2. Work through the **[in-class exercises](Day02_2_InClass_Exercises.md)** (~25 min)
3. Complete the **[homework](Day02_3_Homework.md)** (~45 min)
4. Keep the **[quick reference](Day02_4_Quick_Reference.md)** open while coding
5. Log the day in [`Progress_Log.md`](../../00_Course_Guide/Progress_Log.md)

---

## 📚 What Each File Contains

| File | Purpose | Time | Best For |
|------|---------|------|----------|
| [Day02_1_Lesson](Day02_1_Lesson.md) | Main lesson — **START HERE** | 30 min | Learning the concepts |
| [Day02_2_InClass_Exercises](Day02_2_InClass_Exercises.md) | Guided exercises with collapsible solutions | 25 min | Practice |
| [Day02_3_Homework](Day02_3_Homework.md) | Self-paced homework with solutions | 45 min | Consolidation |
| [Day02_4_Quick_Reference](Day02_4_Quick_Reference.md) | One-page cheat sheet | 5 min | Quick lookup while coding |
| [Day02_5_Extended_Lesson](Day02_5_Extended_Lesson.md) | Long-form version of the lesson with extra examples and built-in exercises | optional | Deeper dive / review |

---

## 🎯 Learning Objectives

By the end of Day 2 you will be able to:

- Use **string interpolation** (`"${sample}.bam"`) and know when single quotes disable it
- Create and manipulate **lists** and **maps** (Python lists and dicts)
- Write **closures** (`{ x -> x * 2 }`, `{ it * 2 }`) — Groovy's lambdas
- Use `collect`, `findAll`, `each` the way you use comprehensions and `filter` in Python
- Recognise Groovy's optional parentheses, which make Nextflow code look "magic"

---

## ✅ Completion Checklist

- [ ] I can interpolate variables with `"${x}"` and know single quotes don't interpolate
- [ ] I know to escape bash variables inside `"""` script blocks (`\$VAR`)
- [ ] I can create and index lists and maps
- [ ] I can write closures with an explicit parameter and with `it`
- [ ] I can translate a Python list comprehension into `collect` / `findAll`
- [ ] I recognise optional parentheses and trailing closures

---

## 🚦 Navigation

⬅️ [Day 1: What Nextflow Actually Is](../Day01_What_Nextflow_Is/README.md) · [📋 Course index](../../README.md) · [Day 3: Your First Nextflow Process](../Day03_First_Process/README.md) ➡️

