# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is a research repository for typhoid immune dynamics modeling. It combines R-based statistical analysis with Python-based documentation site generation using MkDocs Material.

## Key Architecture

- **R-based analysis**: Primary development language for statistical modeling and data analysis
- **Literate programming**: Uses knitr::spin and R Markdown for combining code and documentation
- **Documentation website**: Uses MkDocs Material with automated blog post generation from R scripts
- **Project structure**: 
  - `scratch/` - Exploratory R scripts and prototypes
  - `docs/` - Documentation source files and blog posts
  - `data/` - Raw data files (Excel format)
  - `analysis/` - Structured analysis outputs

## Common Commands

### Documentation Site
```bash
# Install Python dependencies
pip install -r requirements.txt

# Serve documentation site locally
mkdocs serve

# Build documentation site
mkdocs build

# Deploy to GitHub Pages
mkdocs gh-deploy
```

### R Development
This is primarily an RStudio project. Open `typhoid-immune-dynamics.Rproj` in RStudio for proper environment setup. VS Code is also configured in `typhoid-immune-dynamics.code-workspace`. 

Key R helper functions are in `docs/docs_helper_functions.R`:
- `render_blog_post()` - Converts R scripts to blog posts with proper MkDocs YAML headers

Tidyverse code style and package dependencies are preferred, but any tools or style that gets the job done is acceptible.

### Blog Post Workflow
The project uses a custom workflow for converting R analysis scripts into blog posts:
1. Write R scripts with roxygen-style comments (`#'`) for narrative text
    1. To correctly render roxygen-style comments within a function, put two spaces before the comment tag (`  #'`).
2. R scripts should start with a YAML header in roxygen comments:
```r
#' ---
#' title: "Your Title"
#' output:
#'   md_document
#' ---
```
3. Use `render_blog_post()` function to convert to markdown with proper blog metadata:
```r
source('docs/docs_helper_functions.R')
render_blog_post('scratch/your_script.R', categories_list = c('category1', 'category2'))
```
4. Scripts are rendered from `scratch/` to `docs/blog/posts/` automatically
5. It is possible to create folders other than `scratch/` for outputs that are expected to be a more permanent part of the project base. 

## Development Notes

- Uses knitr::spin for literate programming - R scripts with `#'` or `  #'` comments become narrative text
- Blog posts require categories metadata and follow MkDocs Material blog structure
- Site uses custom IDM branding and color scheme
- All R scripts should use 2-space indentation (configured in .Rproj file)
- Git integration for documentation timestamps and contributor attribution