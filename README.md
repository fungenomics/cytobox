# cytobox: README
A toolkit for analyzing single cell RNA-seq data created by the Kleinman Lab. The legacy package,
which is no longer being developed, is available at https://github.com/sjessa/cytobox/

## Installation

You will need to install [`devtools`](https://cran.r-project.org/web/packages/devtools/), and then run:

```r
devtools::install_github("fungenomics/cytobox", build_vignettes = TRUE)

```

## Getting help

Access the barebones vignette from within R for a list of examples:

```r
browseVignettes("cytobox")
```

You can look up the documentation for any `cytobox` function from the R console,
these also contain an example or two for many functions:
```r
?cytobox::tsne
```

Some functions will specify the author in the documentation, whom you could contact directly :)

## Versioning

#### Semantic versioning

`cytobox` uses a semantic versioning scheme, where package versions are of the form x.y.z:

- The patch version z is updated for very minor changes like bug fixes
- The minor version y is udpated for minor but important changes like addition of new functions
- The major version x is updated for major changes, which may introduce backwards incompatibiility

#### Keeping track of cytobox versions in your work

If you use R Markdown, it's a great idea to include the following function as the last chunk in your
documents, which will print a list of all packages loaded, and their versions.

```r
sessionInfo()
```

If you work with scripts, you could print the following to STDOUT/the log files:

```r
packageVersion("cytobox")
```

## Contribution

#### The first time

Clone the repository and create a development branch:
```bash
git clone https://github.com/fungenomics/cytobox.git
cd cytobox
git checkout -b dev-selin
```

Make changes locally on your branch and once the package passes [`R CMD check`](http://r-pkgs.had.co.nz/check.html) with 0 errors,
commit them. Then, push your changes:
```bash
git push -u origin dev-selin
```

If the Travis build passes, there will be a green checkmark next to the 
most recent commit on your branch. When the build has passed, create a pull request on the GitHub website and merge changes in your dev branch with the master branch.

#### From now on

From now on, to continue contributing to the package, switch to your branch and make your
changes there, following the same edit/commit/push workflow. 

If changes have been made to master since first merge the latest changes from the master branch:
```bash
git checkout master
git pull
git checkout dev-selin
# Apply all the changes from master to dev-selin, and set HEAD so that new work
# will happen on top of those changes
git rebase master
```
 
For more on rebasing, see: https://www.atlassian.com/git/tutorials/merging-vs-rebasing
