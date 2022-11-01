# Initialize
rm(list = ls())
wd <- '~/work/gmu_cig/'

startup <- function(path){
  # Loading packages
  library(igraph)
  library(backbone)
  
  # Loading data
  setwd(path)
  dat <- read.csv('cig_student_data.csv')
  courses <- read.csv('cig_courses.csv')
  
  dat$course <- apply(dat[,c('SUBJ', 'NUMB')], 1, paste0, collapse = ' ')
  
  dat <<- dat
  listed <-courses$number
  
  listed <- listed[-grep('/', listed)] # grep('/', listed, value = T) check values and adjust by hand
  listed <<- c(listed, 'EVPP 437', 'BIOL 437', 'EVPP 438', 'BIOL 438', 'EVPP 439', 'BIOL 439', 'GGS 321', 'BIOL 374')
  
  # Functions
  subsetting <<- function(df, program = NULL, year = NULL, listed = NULL){
    if(is.null(program) & is.null(year) & is.null(listed)){
      cat('No subsetting specified. Returning original data.\n')
    } else {
      if(!is.null(program)){df <- df[which(df$PROGR1 == program),]}
      if(!is.null(year)){df <- df[which(floor(df$COURSE_TERM/100) >= year),]}
      if(!is.null(listed)){df <- df[which(df$course %in% listed),]}
    }
    df
  }
  
  data_description <<- function(df){
    length(unique(df$STU_NUMBER)) # unique students: 1247
    
    ncourses <- nrow(unique(df[,c('SUBJ', 'NUMB')])) # 2095 different courses
    
    sort(table(df$SUBJ))          # table of major counts
    plot(density(table(df$SUBJ))) # plotting distribution of courses
  }
  
  toBipartite <<- function(df, subset = 'Median', output = 'igraph'){
    # subsetting to relevant courses by size/count
    course_size <- sort(table(df$course), decreasing = T)
    if(subset == 'Median'){
      topn <- min(which(course_size == round(median(course_size)))) - 1
    } else if(subset == 'All'){
      topn <- length(course_size)
    } else {
      topn <- subset
    }
    topcourses <- course_size[1:topn]
    topnames <- names(topcourses)
    
    df <- df[which(df$course %in% topnames),]
    
    # creating incidence matrix
    mat <- matrix(0, nrow = topn, ncol = length(unique(df$STU_NUMBER)))
    rownames(mat) <- topnames
    colnames(mat) <- unique(df$STU_NUMBER)
    
    for(i in 1:nrow(df)){
      mat[df[i, 'course'], as.character(df[i, 'STU_NUMBER'])] <- mat[df[i, 'course'], as.character(df[i, 'STU_NUMBER'])] + 1
    }
    mat <- ifelse(mat > 0, 1, 0) # don't care about retakes
    
    if(output == 'igraph'){mat <- graph_from_incidence_matrix(mat)}
    mat
  }
  
  commdetect <<- function(bpgraph, resolution_parameter = 1, alpha = 0.05, n_iterations = 2, seed = 1){
    set.seed(seed)
    projCourse <- sdsm(bpg, alpha = 0.1)
    cl <- cluster_leiden(projCourse, resolution_parameter = resolution_parameter, n_iterations = n_iterations)
    
    clusMembership <- membership(cl)
    cat('Inferred', max(clusMembership), 'clusters.', sep = ' ')
    plot(projCourse, vertex.color = rev(kelly.colors(n = 22))[clusMembership])
    list(projCourse, clusMembership)
  }
}


startup(wd)

# Programs: SC-BS-EVSC and LA-BA-EVSS
dsub <- subsetting(dat, program = 'SC-BS-EVSC', year = 2019, listed = listed)
dsub <- subsetting(dat, program = 'LA-BA-EVSS', year = 2018, listed = listed)
dsub <- subsetting(dat, listed = listed)

bpg <- toBipartite(dsub)

comms <- commdetect(bpg, resolution_parameter = 0.2)
