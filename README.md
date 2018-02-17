[![Build Status](https://travis-ci.com/sjessa/cytokit.svg?token=ckZxkx4uN2RZSSwsdpLM&branch=master)](https://travis-ci.com/sjessa/cytokit)

# cytokit
Internal Kleinman Lab toolkit for analyzing single cell RNA-seq data


## Installation

### Hydra/Cedar/Guillimin

To use `cytokit` on Hydra or Guillimin, you can simply load it with `library(cytokit)` in R,
as it's already installed in the Kleinman lab spaces,
and we will take care of keeping the installed version up to date.

### Locally

To install locally, you will need [`devtools`](https://cran.r-project.org/web/packages/devtools/).

Clone the cytokit repository (you will be prompted for your GitHub username and password):
```bash
$ git clone https://github.com/sjessa/cytokit.git
```

And then in R, install the package using devtools:
```
devtools::install_local("path/to/cytokit")
```

To update the package, you will need to pull changes to your copy of the repository:
```bash
$ cd cytokit
$ git pull
```

and run the same `devtools` command in R as above.

<!--
To install locally, you will need to follow a few steps, since the repo is private:

1. Install [`devtools`](https://cran.r-project.org/web/packages/devtools/)
2. Generate a personal access token on Github here: https://github.com/settings/tokens  
    Generate new token > Type in "cytokit" in *Token description* > Tick the checkbox for "repo" > Generate token (green button)
3. Copy the token (it will be a string of letters and numbers, e.g. `abc123`)
4. Install `cytokit` and include the token as an argument:

```r
devtools::install_github("sjessa/cytokit", auth_token = "abc123")

``` -->
