# 🧠 Linux Basics

This section covers the following topics:
1. Why we are using Linux  
2. How to install Ubuntu  
3. Basic commands and their importance  
4. Special tricks  
5. How to create an environment  

---

## 🧩 1️⃣ Why We Are Using Linux
- Linux can handle large data efficiently, which is essential in bioinformatics and NGS workflows.  
- It is stable, secure, and supports a wide range of command-line tools used in genomics.  
- Two popular Linux distributions (subsystems) are **Ubuntu** and **Kali**.  
- Linux organizes all directories in a **tree-like structure**.  
- In Linux, **folders are called directories**.  

---

## 💻 2️⃣ How to Install Ubuntu
I downloaded **Ubuntu 22.04.5 LTS** from the Microsoft Store.  
While setting it up, I created a **username (ID)** and a **password**.  

> ⚠️ The password is very important. If forgotten, Ubuntu must be uninstalled and reinstalled.  

After installation, the commands were practiced directly in Ubuntu.

---

## 📘 3️⃣ Basic Commands and Their Importance

| Command | Description |
|----------|-------------|
| `ls` | List directory contents |
| `grep` | Print lines that match a pattern |
| `mv` | Move or rename files |
| `cp` | Copy files and directories |
| `rm` | Remove a file (⚠️ Never use `rm -r` casually; see note below) |
| `cat` | Concatenate files and print on standard output |
| `more` | View file contents one screen at a time |
| `less` | View file contents (use `q` to quit) |
| `wc` | Print newline, word, and byte counts for each file |
| `head` | Display the first part of files |
| `tail` | Display the last part of files |
| `clear` | Clear the terminal screen |
| `pwd` | Print the current working directory |
| `realpath` | Print the full resolved path |

---

### ⚠️ Important Note on `rm -r`
- The command `rm -r` means **remove recursively** — it deletes a folder and everything inside it (subfolders and files).  
- It is **permanent** and cannot be undone.  
- Always double-check before using it.  


✅ **Safe usage examples:**  


# Remove a single file
rm filename

# Remove a folder interactively (asks for confirmation)
rm -ri foldername

## ⚙️ 4️⃣  Special Tricks
🔹 Navigating Directories
|--------|------------|
| `cd` | change directory (Always put a space after cd) |
| `/` | represents the root directory |
| `ls` | lists contents |
| `cd ..` | move one directory up |
| `cd ../..` | move two directories up |
| `clear` | clear all written commands and start fresh |
| `pwd` | shows your full working path |
| `realpath` | shows the absolute path of a file or folder |

Pressing Ctrl + Z stops the currently running command and brings you back to the user prompt.

## 🔹 Tab Completion
After typing part of a directory or file name, press Tab to auto-complete it.
Example: typing `cd Int` + pressing Tab completes to
`cd Introduction_to_Unix_Shell_practice`

## 🔹 Using Help and Flags
Flags modify how commands behave.

ls -a     # show hidden files
ls -t     # sort by time modified
ls -lh    # list in human-readable format with permissions
ls -lhtr  # detailed listing, sorted by time (reverse order)

🔹 Creating and Managing Files/Folders
| Action | Command |
|--------|----------|
| Make a new directory | `mkdir newdirectory` | 
| Move a file to a folder |	`mv file newfolder/` |
| Copy a file to a folder |	`cp file newfolder/` |
| Rename a file |	`mv oldname newname` |
 Copy and rename a file |	`cp oldname newname` |

🔹 Viewing Text Files
| Action | Command |
|---------|---------|
| View full file |	`less filename` (press q to quit) |
| View first 10 lines |	`head filename` |
| View last 10 lines |	`tail filename` |
| View first N lines |	`head -n [N] filename` |
| View last N lines |	`tail -n [N] filename` |
| View a particular line |	`head -n [N] filename | tail -n 1` |

🔹 Using cat and Redirection
Purpose	Command
Display file(s)	cat file1 file2
Combine multiple files into a new file	cat file1 file2 > newfile
Append to an existing file	cat file3 >> newfile

### ⚠️ Be careful:

Using `>` overwrites existing content.

Using `>>` adds content to the existing file.

Empty files cannot be appended using `>>`.

🔹 Using grep
To print lines that contain a specific word or pattern:

```bash
grep ATOM filename.pdb
This shows only those lines where “ATOM” appears.

🔹 Editing Files Using nano
nano is a simple text editor used directly in the terminal.

| Action | Command |
| ------ | -------- |
| Open or create a file |	`nano filename` |
| Save changes |	Press Ctrl + O, then Enter |
| Exit | nano	Press Ctrl + X |

Example:

```bash
nano notes.txt
This opens a text file called notes.txt. You can write, edit, then save and close it.

## 🌿 5️⃣ Environment Setup
🔹 What Is an Environment?
- An environment is a separate workspace in which specific tools or packages are installed.
- I created different environments for different tools because:
1. It helps to locate tools easily.
2. Some tools cannot run in the same environment due to dependency conflicts.

🔹 Commands for Managing Environments
| Action | Command |
| Create a new environment | `mamba create -n envname` |
| Activate an environment |	`mamba activate envname` |
| Deactivate the current environment |	`mamba deactivate` |
| List all environments |	`mamba env list` |
| To view tools installed in an environment | `mamba list` |
| View packages in an environment |	`mamba list -n envname` |
| Remove an environment |	`mamba env remove -n envname` | 

