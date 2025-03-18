# docs_helper_functions.R
# functions to help render the docs automatically from R code

render_blog_post = function(input_file,
                            categories_list,
                            date_created=lubridate::today(),
                            authors_list='famulare',
                            draft = FALSE,
                            knit_root_dir = rprojroot::find_rstudio_root_file(),
                            debug_blog_header = FALSE){
  # this function renders any commented .R or .Rmd files into a blog post, with
  # files put in the right location in the docs: /docs/blog/posts
  #
  # required arguments
  # - input_file:  input file name with path relative to working directory
  # - categories_list: list of categories for the blog post. Categories_list is required and has no default. 
  #                    Please include at least one valid category for the blog post.
  # 
  # optional arguments
  # - date_created = today(): date post is created. default is today().
  # - authors_list = 'famulare': list of github handle(s) of post author(s).
  # - draft = FALSE: TRUE/FALSE for if post is draft status or not
  # - knit_root_dir = rprojroot::find_rstudio_root_file(): root directory from which relative paths are defined
  # - output_dir = './docs/blog/posts': default directory relative to rprojroot
  # - debug_blog_header = FALSE: if TRUE, output the constructed blog header and do not render.
  # 
  # OUTPUT
  # - if debug_blog_header == FALSE: then blog post rendered for mkdocs and put in correct directory
  # - else if debug_blog_header == TRUE: return blog post header string
  #
  
  # input handling
    # categories_list is required and cannot be empty
    if (length(categories_list)==0){
      simpleError('categories_list is required and has no default. Please include at least one valid category for the blog post.')
    }
  
    date_created = as.character(date_created)
  
  blog_header = paste(sep='\n',
    '---',
    paste0('draft: ',tolower(as.character(draft))),
    'date:',
    paste0('  created: ', date_created),
    'authors:',
    paste('  - ',authors_list,collapse='\n', sep=''),
    'categories:',
    paste('  - ',categories_list,collapse='\n', sep=''),
    '---\n'
    )

  if (debug_blog_header){
    return(blog_header)
  } else {
    
    # get output .md to modify for mkdocs material blog specification
    output_dir = paste0(knit_root_dir,'/docs/blog/posts')
    output_file <- paste(output_dir,'/',sub("^(?:.*/)?([^/]+)\\.[^.]+$", "\\1", input_file),'.md',sep='')
    
    # render knitr::spin commented or Rmd file
    rmarkdown::render(input_file,
                      knit_root_dir = knit_root_dir,
                      output_dir = output_dir)
    
    # read blog post
    blog_post = readLines(output_file)

    # fixing the fact the render uses absolute paths when the output directory isn't the input directory
    # https://stackoverflow.com/questions/70098862/force-relative-paths-in-knitrinclude-graphics
    blog_post = gsub(paste0("(?<=!\\[]\\()[^)]*?(?=",paste0(paste0(sub("^(?:.*/)?([^/]+)\\.[^.]+$", "\\1", input_file),'_files')),")"), "", blog_post, perl = TRUE)
    
    # append MkDocs material blog YAML header
    writeLines(c(blog_header,blog_post),output_file)
  }
}


