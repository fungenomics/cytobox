[![Build Status](https://travis-ci.com/sjessa/cytokit.svg?token=ckZxkx4uN2RZSSwsdpLM&branch=master)](https://travis-ci.com/sjessa/cytokit)

# cytokit
Internal Kleinman Lab toolkit for analyzing single cell RNA-seq data


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

## Getting help

Access the very barebones vignette from within R for a list of examplesx:

```r
browseVignettes("cytokit")
```

You can look up the documentation for any `cytokit` function from the R console:
```r
?cytokit::ggColours
```

Some functions will specify the author in the documentation, whom you could contact directly :)

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

### From now on

From now on, to continue contributing to the package, switch to your branch and make your
changes there, following the same edit/commit/push workflow. 

If changes have been made to master since first merge the latest changes from the master branch:
```bash
git checkout master
git pull
git checkout dev-selin
git merge master # Merge changes from master branch onto dev-selin
```


(There is a more sophisticated way of doing this merge, called rebasing: https://www.atlassian.com/git/tutorials/merging-vs-rebasing)
