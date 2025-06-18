---
layout: page
icon: fas fa-stream
order: 1
---


> This guide will walk you through the process of installing R and the BIGr package on your computer.
{: .prompt-info }

---

## Installing R

R is a powerful, open-source programming language and software environment for statistical computing and graphics. You will need to install R before you can use the BIGr package.

Please follow the instructions for your operating system.

### Windows

1.  **Download R:** Go to the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/bin/windows/base/) website.
2.  **Run the Installer:** Click on the "Download R-X.X.X for Windows" link (the X's will be the latest version number). Run the downloaded `.exe` file.
3.  **Follow the Setup Wizard:**
    * Select your preferred language.
    * Read the license agreement and click "Next".
    * Choose an installation location or accept the default.
    * Select the components to install (the defaults are usually sufficient).
    * Accept the startup options or customize them if you are an advanced user.
    * Choose a Start Menu folder.
    * Select any additional tasks, like creating a desktop icon.
    * Click "Install" and wait for the installation to complete.
    * Click "Finish" to exit the setup wizard.

### macOS

1.  **Download R:** Go to the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/bin/macosx/) website.
2.  **Choose the Correct Package:**
    * For Macs with Apple Silicon (M1, M2, etc.), download the package for "Apple Silicon (arm64)".
    * For older Macs with Intel processors, download the package for "Intel 64-bit".
3.  **Run the Installer:** Open the downloaded `.pkg` file.
4.  **Follow the Installation Steps:**
    * Click "Continue" through the welcome and license agreement screens.
    * Agree to the software license agreement.
    * Select the installation destination and click "Install". You may be prompted to enter your password.
    * Once the installation is complete, you can close the installer.

### Linux

The installation process for Linux can vary depending on your distribution.

### Verify Installation

To verify that R has been installed correctly, you can open a terminal (or the R console on Windows and macOS) and type:
```R
R --version
```
This should print the installed R version.

## Installing the BIGr R Package
Once you have R installed, you can install the BIGr package directly from the R console.

There are two primary methods for installing the BIGr package: from CRAN or from GitHub for the latest development version.

> You will likely be prompted in the terminal to install the package dependencies or update any existing packages. You must select one of the options before the installation will proceed.
{: .prompt-tip }

### Option 1: Install from CRAN (Comprehensive R Archive Network)
This is the recommended method for most users as it provides the latest stable version of the package.

1. Open your R console.
2. Run the following command:
```R
install.packages("BIGr")
```

### Option 2: Install the Development Version from GitHub
If you need the absolute latest features and are comfortable with potentially less stable code, you can install the package from the Breeding Insight GitHub repository.

1. Open your R console.
2. You will first need to install the `remotes`  and `BiocManager` packages if you don't already have them:
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    install.packages("remotes")
```
3. Install the BIGr package from GitHub:
```R
BiocManager::install("Breeding-Insight/BIGr", dependencies = TRUE)
```

After the installation is complete, you can load the BIGr package into your R session to start using it:
```R
library(BIGr)
```

> You are now ready to use the functionalities provided by the BIGr package!
{: .prompt-info }