[![Build Status](https://travis-ci.com/sjessa/cytokit.svg?token=ckZxkx4uN2RZSSwsdpLM&branch=master)](https://travis-ci.com/sjessa/cytokit)

# cytokit: README
Kleinman Lab toolkit for analyzing single cell RNA-seq data.
<!---Link to vignette: [https://rawgit.com/sjessa/cytokit/master/vignettes/cytokit.html](https://rawgit.com/sjessa/cytokit/master/vignettes/cytokit.html) -->


## Installation

You will need to install [`devtools`](https://cran.r-project.org/web/packages/devtools/), and then run:

```r
devtools::install_github("sjessa/cytokit", build_vignettes = TRUE)

```

## Getting help

Access the barebones vignette from within R for a list of examples:

```r
browseVignettes("cytokit")
```

<!---Or checkout the version saved in the repository (not necessarily up to date!) here: https://rawgit.com/sjessa/cytokit/master/vignettes/cytokit.html-->

You can look up the documentation for any `cytokit` function from the R console,
these also contain an example or two for many functions:
```r
?cytokit::tsne
```

Some functions will specify the author in the documentation, whom you could contact directly :)

## Versioning

#### Semantic versioning

`cytokit` uses a semantic versioning scheme, where package versions are of the form x.y.z:

- The patch version z is updated for very minor changes like bug fixes
- The minor version y is udpated for minor but important changes like addition of new functions
- The major version x is updated for major changes, which may introduce backwards incompatibiility

#### Keeping track of cytokit versions in your work

If you use R Markdown, it's a great idea to include the following function as the last chunk in your
documents, which will print a list of all packages loaded, and their versions.

```r
sessionInfo()
```

If you work with scripts, you could print the following to STDOUT/the log files:

```r
packageVersion("cytokit")
```

## Contribution

#### The first time

Clone the repository and create a development branch:
```bash
git clone https://github.com/sjessa/cytokit.git
cd cytokit
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
