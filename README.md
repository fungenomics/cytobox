[![Build Status](https://travis-ci.com/sjessa/cytokit.svg?token=ckZxkx4uN2RZSSwsdpLM&branch=master)](https://travis-ci.com/sjessa/cytokit)

# cytokit: README
Internal Kleinman Lab toolkit for analyzing single cell RNA-seq data.
Link to vignette: [https://rawgit.com/sjessa/cytokit/dev-selin/vignettes/cytokit.html](https://rawgit.com/sjessa/cytokit/dev-selin/vignettes/cytokit.html)


## Installation

#### Hydra/Cedar/Guillimin

To use `cytokit` on Hydra or Guillimin, you can simply load it with `library(cytokit)` in R,
as it's already installed in the Kleinman lab spaces. To update the version on the server, you can use the instructions below.

#### Locally

Since the repo is private, there are a couple steps to install `cytokit`:

1. Install [`devtools`](https://cran.r-project.org/web/packages/devtools/)
2. Generate a personal access token on Github here: https://github.com/settings/tokens  
    Generate new token > Type in "cytokit" in *Token description* > Tick the checkbox for "repo" > Generate token (green button)  
   Your token will be a string of letters and numbers, e.g. `abc123`. (For more details on these tokens, see the [GitHub documentation](https://help.github.com/articles/creating-a-personal-access-token-for-the-command-line/))
3. Install `cytokit` and include the token as an argument:

```r
devtools::install_github("sjessa/cytokit", auth_token = "abc123", build_vignettes = TRUE)

```

You'll need this token again in the future to update `cytokit`, so it's a good
idea to save the above command in a script, or save your token somewhere you 
can copy-paste it. 

#### Updating

To update `cytokit`, re-run the install command in Step 3 above.

If you don't know your old token, go to https://github.com/settings/tokens, click 'Edit' on
your cytokit token, click 'Regenerate' and re-run the command with the new token.

## Getting help

Access the barebones vignette from within R for a list of examples:

```r
browseVignettes("cytokit")
```

Or checkout the version saved in the repository (not necessarily up to date!) here: https://rawgit.com/sjessa/cytokit/dev-selin/vignettes/cytokit.html

You can look up the documentation for any `cytokit` function from the R console,
these also contain an example or two for many functions:
```r
?cytokit::tsne
```

Some functions will specify the author in the documentation, whom you could contact directly :)

## Versioning

#### Semantic versioning

`cytokit` uses a semantic versioning scheme, where package versions are of the form x.y.z:

- The z version is updated for very minor changes like bug fixes
- The minor version y is udpated for minor but important changes like addition of new functions
- The major version z is updated for major changes, which may introduce backwards incompatibiility

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
