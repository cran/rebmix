setGeneric("RNGMIX",
  function(model = "RNGMIX",
    Dataset.name = character(),
    rseed = -1,
    n = numeric(),
    Theta = list(), ...)
  standardGeneric("RNGMIX"))

setGeneric("REBMIX",
  function(model = "REBMIX",
    Dataset = list(),
    Preprocessing = character(),
    cmax = 15,
    cmin = 1,
    Criterion = "AIC",
    pdf = character(),
    theta1 = numeric(),
    theta2 = numeric(),
    theta3 = numeric(),
    K = "auto",
    ymin = numeric(),
    ymax = numeric(),
    ar = 0.1,
    Restraints = "loose",
    Mode = "outliersplus",
### Panic Branislav.    
    EMcontrol = NULL, ...)
### End    
  standardGeneric("REBMIX"))

### Panic Branislav.
  setGeneric("EMMIX",
  function(model = "REBMIX",
    Dataset = list(),
    Theta = NULL,
    EMcontrol = NULL, ...)
  standardGeneric("EMMIX"))
### End

setGeneric(".IC", function(x = NULL, Criterion = "AIC", pos = 1, ...) standardGeneric(".IC"))

setGeneric("logL", function(x = NULL, pos = 1, ...) standardGeneric("logL"))
setGeneric("AIC", function(x = NULL, pos = 1, ...) standardGeneric("AIC"))
setGeneric("AIC3", function(x = NULL, pos = 1, ...) standardGeneric("AIC3"))
setGeneric("AIC4", function(x = NULL, pos = 1, ...) standardGeneric("AIC4"))
setGeneric("AICc", function(x = NULL, pos = 1, ...) standardGeneric("AICc"))
setGeneric("BIC", function(x = NULL, pos = 1, ...) standardGeneric("BIC"))
setGeneric("CAIC", function(x = NULL, pos = 1, ...) standardGeneric("CAIC"))
setGeneric("HQC", function(x = NULL, pos = 1, ...) standardGeneric("HQC"))
setGeneric("MDL2", function(x = NULL, pos = 1, ...) standardGeneric("MDL2"))
setGeneric("MDL5", function(x = NULL, pos = 1, ...) standardGeneric("MDL5"))
setGeneric("AWE", function(x = NULL, pos = 1, ...) standardGeneric("AWE"))
setGeneric("CLC", function(x = NULL, pos = 1, ...) standardGeneric("CLC"))
setGeneric("ICL", function(x = NULL, pos = 1, ...) standardGeneric("ICL"))
setGeneric("ICLBIC", function(x = NULL, pos = 1, ...) standardGeneric("ICLBIC"))
setGeneric("PRD", function(x = NULL, pos = 1, ...) standardGeneric("PRD"))
setGeneric("SSE", function(x = NULL, pos = 1, ...) standardGeneric("SSE"))
setGeneric("PC", function(x = NULL, pos = 1, ...) standardGeneric("PC"))

setGeneric("demix",
  function(x = NULL,
    pos = 1,
    variables = expression(1:d), ...)
  standardGeneric("demix"))

setGeneric("pemix",
  function(x = NULL,
    pos = 1,
    variables = expression(1:d),
    lower.tail = TRUE,
    log.p = FALSE, ...)
  standardGeneric("pemix"))

setGeneric("dfmix",
  function(x = NULL,
    Dataset = data.frame(),
    pos = 1,
    variables = expression(1:d), ...)
  standardGeneric("dfmix"))

setGeneric("pfmix",
  function(x = NULL,
    Dataset = data.frame(),
    pos = 1,
    variables = expression(1:d),
    lower.tail = TRUE,
    log.p = FALSE, ...)
  standardGeneric("pfmix"))

setGeneric("split",
  function(p = NULL,
    Dataset = data.frame(),
    class = numeric(), ...)
  standardGeneric("split"))

setGeneric("chunk",
  function(x = NULL,
    variables = expression(1:d))
  standardGeneric("chunk"))

setGeneric("boot",
  function(x = NULL,
    rseed = -1,
    pos = 1,
    Bootstrap = "parametric",
    B = 100,
    n = numeric(),
    replace = TRUE,
    prob = numeric(), ...)
  standardGeneric("boot"))

setGeneric("RCLRMIX",
  function(model = "RCLRMIX",
    x = NULL,
    Dataset = NULL,
    pos = 1,
    Zt = factor(), 
    Rule = character(), ...)
  standardGeneric("RCLRMIX"))

setGeneric("RCLSMIX",
  function(model = "RCLSMIX",
    x = list(),
    Dataset = data.frame(),
    Zt = factor(), ...)
  standardGeneric("RCLSMIX"))

setGeneric("BFSMIX",
  function(model = "RCLSMIX",
    x = list(),
    Dataset = data.frame(),
    Zt = factor(), ...)
  standardGeneric("BFSMIX"))

setGeneric("a.d", function(x = NULL) standardGeneric("a.d"))
setGeneric("a.theta1<-", function(x = NULL, l = numeric(), value = numeric()) standardGeneric("a.theta1<-"))
setGeneric("a.theta2<-", function(x = NULL, l = numeric(), value = numeric()) standardGeneric("a.theta2<-"))
setGeneric("a.theta3<-", function(x = NULL, l = numeric(), value = numeric()) standardGeneric("a.theta3<-"))
### Panic Branislav.
setGeneric("a.theta1.all<-", function(x = NULL, value = numeric()) standardGeneric("a.theta1.all<-"))
setGeneric("a.theta2.all<-", function(x = NULL, value = numeric()) standardGeneric("a.theta2.all<-"))
setGeneric("a.theta3.all<-", function(x = NULL, value = numeric()) standardGeneric("a.theta3.all<-"))
### End
setGeneric("a.Dataset.name", function(x = NULL) standardGeneric("a.Dataset.name"))
setGeneric("a.rseed", function(x = NULL) standardGeneric("a.rseed"))
setGeneric("a.n", function(x = NULL) standardGeneric("a.n"))
setGeneric("a.Theta", function(x = NULL, pos = 0) standardGeneric("a.Theta"))
setGeneric("a.Dataset", function(x = NULL, pos = 0, ...) standardGeneric("a.Dataset"))
setGeneric("a.Zt", function(x = NULL) standardGeneric("a.Zt"))
setGeneric("a.w", function(x = NULL, pos = 0) standardGeneric("a.w"))
### Panic Branislav.
setGeneric("a.w<-", function(x = NULL, value = numeric()) standardGeneric("a.w<-"))
### End
setGeneric("a.Variables", function(x = NULL) standardGeneric("a.Variables"))
setGeneric("a.ymin", function(x = NULL) standardGeneric("a.ymin"))
setGeneric("a.ymax", function(x = NULL) standardGeneric("a.ymax"))

setGeneric("a.Preprocessing", function(x = NULL) standardGeneric("a.Preprocessing"))
setGeneric("a.cmax", function(x = NULL) standardGeneric("a.cmax"))
setGeneric("a.cmin", function(x = NULL) standardGeneric("a.cmin"))
setGeneric("a.Criterion", function(x = NULL) standardGeneric("a.Criterion"))
setGeneric("a.pdf", function(x = NULL) standardGeneric("a.pdf"))
setGeneric("a.theta1", function(x = NULL) standardGeneric("a.theta1"))
setGeneric("a.theta2", function(x = NULL) standardGeneric("a.theta2"))
setGeneric("a.theta3", function(x = NULL) standardGeneric("a.theta3"))
setGeneric("a.theta1.all", function(x = NULL, pos = 1) standardGeneric("a.theta1.all"))
setGeneric("a.theta2.all", function(x = NULL, pos = 1) standardGeneric("a.theta2.all"))
setGeneric("a.theta3.all", function(x = NULL, pos = 1) standardGeneric("a.theta3.all"))
setGeneric("a.K", function(x = NULL) standardGeneric("a.K"))
setGeneric("a.K<-", function(x = NULL, value = numeric()) standardGeneric("a.K<-"))
setGeneric("a.y0", function(x = NULL) standardGeneric("a.y0"))
setGeneric("a.ar", function(x = NULL) standardGeneric("a.ar"))
setGeneric("a.Restraints", function(x = NULL) standardGeneric("a.Restraints"))
setGeneric("a.Mode", function(x = NULL) standardGeneric("a.Mode"))
setGeneric("a.summary", function(x = NULL, col.name = character(), pos = 0) standardGeneric("a.summary"))
setGeneric("a.pos", function(x = NULL) standardGeneric("a.pos"))
setGeneric("a.opt.c", function(x = NULL) standardGeneric("a.opt.c"))
setGeneric("a.opt.IC", function(x = NULL) standardGeneric("a.opt.IC"))
setGeneric("a.opt.logL", function(x = NULL) standardGeneric("a.opt.logL"))
setGeneric("a.opt.Dmin", function(x = NULL) standardGeneric("a.opt.Dmin"))
setGeneric("a.opt.D", function(x = NULL) standardGeneric("a.opt.D"))
setGeneric("a.all.K", function(x = NULL) standardGeneric("a.all.K"))
setGeneric("a.all.IC", function(x = NULL) standardGeneric("a.all.IC"))

setGeneric("a.Bootstrap", function(x = NULL) standardGeneric("a.Bootstrap"))
setGeneric("a.B", function(x = NULL) standardGeneric("a.B"))
setGeneric("a.replace", function(x = NULL) standardGeneric("a.replace"))
setGeneric("a.prob", function(x = NULL) standardGeneric("a.prob"))
setGeneric("a.c", function(x = NULL) standardGeneric("a.c"))
setGeneric("a.c.se", function(x = NULL) standardGeneric("a.c.se"))
setGeneric("a.c.cv", function(x = NULL) standardGeneric("a.c.cv"))
setGeneric("a.c.mode", function(x = NULL) standardGeneric("a.c.mode"))
setGeneric("a.c.prob", function(x = NULL) standardGeneric("a.c.prob"))
setGeneric("a.w.se", function(x = NULL) standardGeneric("a.w.se"))
setGeneric("a.w.cv", function(x = NULL) standardGeneric("a.w.cv"))
setGeneric("a.Theta.se", function(x = NULL) standardGeneric("a.Theta.se"))
setGeneric("a.Theta.cv", function(x = NULL) standardGeneric("a.Theta.cv"))

setGeneric("a.Zp", function(x = NULL, s = expression(c)) standardGeneric("a.Zp"))
setGeneric("a.p", function(x = NULL, s = expression(c)) standardGeneric("a.p"))
setGeneric("a.pi", function(x = NULL, s = expression(c)) standardGeneric("a.pi"))
setGeneric("a.P", function(x = NULL, s = expression(c)) standardGeneric("a.P"))
setGeneric("a.tau", function(x = NULL, s = expression(c)) standardGeneric("a.tau"))
setGeneric("a.from", function(x = NULL) standardGeneric("a.from"))
setGeneric("a.to", function(x = NULL) standardGeneric("a.to"))
setGeneric("a.EN", function(x = NULL) standardGeneric("a.EN"))
setGeneric("a.ED", function(x = NULL) standardGeneric("a.ED"))
setGeneric("a.Rule", function(x = NULL) standardGeneric("a.Rule"))

setGeneric("a.o", function(x = NULL) standardGeneric("a.o"))
setGeneric("a.s", function(x = NULL) standardGeneric("a.s"))
setGeneric("a.ntrain", function(x = NULL) standardGeneric("a.ntrain"))
setGeneric("a.Zr", function(x = NULL) standardGeneric("a.Zr"))
setGeneric("a.ntest", function(x = NULL) standardGeneric("a.ntest"))
setGeneric("a.CM", function(x = NULL) standardGeneric("a.CM"))
setGeneric("a.Accuracy", function(x = NULL) standardGeneric("a.Accuracy"))
setGeneric("a.Error", function(x = NULL) standardGeneric("a.Error"))
setGeneric("a.Precision", function(x = NULL) standardGeneric("a.Precision"))
setGeneric("a.Sensitivity", function(x = NULL) standardGeneric("a.Sensitivity"))
setGeneric("a.Specificity", function(x = NULL) standardGeneric("a.Specificity"))
setGeneric("a.Chunks", function(x = NULL) standardGeneric("a.Chunks"))

setGeneric("a.levels", function(x = NULL) standardGeneric("a.levels"))
setGeneric("a.train", function(x = NULL) standardGeneric("a.train"))
setGeneric("a.test", function(x = NULL) standardGeneric("a.test"))

### Panic Branislav.
setGeneric("a.strategy", function(x = NULL) standardGeneric("a.strategy"))
setGeneric("a.strategy<-", function(x = NULL, value = character()) standardGeneric("a.strategy<-"))
setGeneric("a.variant", function(x = NULL) standardGeneric("a.variant"))
setGeneric("a.variant<-", function(x = NULL, value = character()) standardGeneric("a.variant<-"))
setGeneric("a.acceleration", function(x = NULL) standardGeneric("a.acceleration"))
setGeneric("a.acceleration<-", function(x = NULL, value = character()) standardGeneric("a.acceleration<-"))
setGeneric("a.tolerance", function(x = NULL) standardGeneric("a.tolerance"))
setGeneric("a.tolerance<-", function(x = NULL, value = numeric()) standardGeneric("a.tolerance<-"))
setGeneric("a.acceleration.multiplier", function(x = NULL) standardGeneric("a.acceleration.multiplier"))
setGeneric("a.acceleration.multiplier<-", function(x = NULL, value = numeric()) standardGeneric("a.acceleration.multiplier<-"))
setGeneric("a.maximum.iterations", function(x = NULL) standardGeneric("a.maximum.iterations"))
setGeneric("a.maximum.iterations<-", function(x = NULL, value = numeric()) standardGeneric("a.maximum.iterations<-"))
setGeneric("a.summary.EM", function(x = NULL, col.name = character(), pos = 0) standardGeneric("a.summary.EM"))
setGeneric("a.eliminate.zero.components", function(x = NULL) standardGeneric("a.eliminate.zero.components"))
setGeneric("a.eliminate.zero.components<-", function(x = NULL, value = numeric()) standardGeneric("a.eliminate.zero.components<-"))
### End

setGeneric("a.Y", function(x = NULL) standardGeneric("a.Y"))
setGeneric("a.Y<-", function(x = NULL, value = numeric()) standardGeneric("a.Y<-"))
setGeneric("a.h", function(x = NULL) standardGeneric("a.h"))
setGeneric("a.ns", function(x = NULL) standardGeneric("a.ns"))

### Panic Branislav & Marko Nagode.
setGeneric("optbins",
  function(Dataset = list(),
    Rule = "Knuth equal",  
    ymin = numeric(),
    ymax = numeric(), 
    kmin = numeric(),
    kmax = numeric(), ...)
  standardGeneric("optbins"))
  
setGeneric("bins",
  function(Dataset = list(),
    K = matrix(),
    ymin = numeric(),
    ymax = numeric(), ...)
  standardGeneric("bins"))  
### End

setGeneric("fhistogram",
  function(x = NULL,
    Dataset = data.frame(),
    K = numeric(),
    ymin = numeric(),
    ymax = numeric(), 
    shrink = FALSE, ...)
  standardGeneric("fhistogram"))
  
setGeneric("chistogram",
  function(x = NULL,
    Dataset = data.frame(),
    K = numeric(),
    ymin = numeric(),
    ymax = numeric(), ...)
  standardGeneric("chistogram"))
  
setGeneric("mapclusters",
  function(x = NULL,
    Dataset = data.frame(),
    s = expression(c), ...)
  standardGeneric("mapclusters"))

setGeneric("labelmoments",
  function(Zp = array(),
    cmax = integer(), 
    Sigma = 1.0, ...)
  standardGeneric("labelmoments"))
  
setGeneric("mergelabels",
  function(A = list(),
    w = numeric(), 
    k = 2, ...)
  standardGeneric("mergelabels"))  

