# Workspace Index

> Records of all AI Agent work records across all developers

---

## Overview

This directory tracks records for all developers working with AI Agents on this project.

### File Structure

```
workspace/
|-- index.md              # This file - main index
+-- {developer}/          # Per-developer directory
    |-- index.md          # Personal index with session history
    |-- journal-N.md     # Journal files (sequential: 1, 2, 3...)
```

---

## Active Developers

| Developer | Last Active | Sessions | Active File |
|-----------|-------------|----------|-------------|
| xhhuang24 | 2026-05-27 | 1 | journal-1.md |

---

## Getting Started

### For New Developers

Run the initialization script:

```bash
python3 ./.trellis/scripts/init_developer.py <your-name>
```

This will:
1. Create your identity file (gitignored)
2. Create your progress directory
3. Create your personal index
4. Create initial journal file

### For Returning Developers

1. Get your developer name:
   ```bash
   python3 ./.trellis/scripts/get_developer.py
   ```

2. Read your personal index:
   ```bash
   cat .trellis/workspace/$(python3 ./.trellis/scripts/get_developer.py)/index.md
   ```

---

## Guidelines

### Journal File Rules

- **Max 2000 lines** per journal file
- When limit is reached, create `journal-{N+1}.md`
- Update your personal `index.md` when creating new files

---

**Language**: All documentation must be written in **English**.
