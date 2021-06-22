classroom <- read.csv("http://www-personal.umich.edu/~bwest/classroom.csv")

schoolLookupVec <- unique(classroom[c("classid","schoolid")])[,"schoolid"]

dat <- with(classroom,
            list(Ni           = length(unique(childid)),
                 Nj           = length(unique(classid)),
                 Nk           = length(unique(schoolid)),
                 classid      = classid,
                 schoolid     = schoolid,
                 schoolLookup = schoolLookupVec,
                 mathgain     = mathgain))

resStan <- stan(model_code = stan_code, data = dat, chains = 4, iter = 10000, warmup = 1000, thin = 10)