#====================================================================
# 20260429EJ
# MSE template for a full feedback model with a4a sca
# Multistock MP
#====================================================================

#--------------------------------------------------------------------
# load libraries and functions
#--------------------------------------------------------------------
library(FLa4a)
library(FLBRP)
library(ggplotFL)
library(mse)
library(msemodules)
load("setup.RData")
source("utilities.R")

#--------------------------------------------------------------------
# load OMs, OEMs and IEMs
#--------------------------------------------------------------------
load("../data/hke811_conditioned.RData")
load("../data/dps811_conditioned.RData")
load("../data/mut9_conditioned.RData")
load("../data/nep9_conditioned.RData")
load("../data/nep11_conditioned.RData")
plot(window(FLStocks(hke811=hke811.om@stock, dps811=dps811.om@stock, mut9=mut9.om@stock,
                     nep9=nep9.om@stock, nep11=nep11.om@stock), end=2034))

#====================================================================
# MP
#====================================================================
hke811.rule <- mpCtrl(
  # estimation method: full feedback (a4a) sca
  est = mseCtrl(method=hke811.sa, args=list(fmodel=fmodel(hke811.fit),
                                            qmodel=lapply(qmodel(hke811.fit),"formula"),
                                            srmodel=srmodel(hke811.fit), dfm=hke811.dfm, update=TRUE)),
  # parametrizing the HCR
  phcr = mseCtrl(method=wmm_ca4.phcr, args=list(interval=rpyi, fmsy_proxy="fmax", blim_proxy="fmax", blim_proxy_multiplier = 0.25)),
  # hcr: westmed MAP
  hcr = mseCtrl(method=westmedmap.hcr),
  # (i)mplementation (sys)tem: effort (in fact F if relationship E~F is 1)
  isys = mseCtrl(method=f.is)
)

dps811.rule <- mpCtrl(
  # estimation method: full feedback (a4a) sca
  est = mseCtrl(method=dps811.sa, args=list(fmodel=fmodel(dps811.fit),
                                            qmodel=lapply(qmodel(dps811.fit),"formula"),
                                            srmodel=srmodel(dps811.fit), dfm=dps811.dfm, update=TRUE)),
  # parametrizing the HCR
  phcr = mseCtrl(method=wmm_ca4.phcr, args=list(interval=rpyi, fmsy_proxy="f0.1", blim_proxy="f0.1", blim_proxy_multiplier = 0.25)),
  # hcr: westmed MAP
  hcr = mseCtrl(method=westmedmap.hcr),
  # (i)mplementation (sys)tem: effort (in fact F if relationship E~F is 1)
  isys = mseCtrl(method=f.is)
)

mut9.rule <- mpCtrl(
  # estimation method: full feedback (a4a) sca
  est = mseCtrl(method=mut9.sa, args=list(fmodel=fmodel(mut9.fit),
                                          qmodel=lapply(qmodel(mut9.fit),"formula"),
                                          srmodel=srmodel(mut9.fit), dfm=mut9.dfm, update=TRUE)),
  # parametrizing the HCR
  phcr = mseCtrl(method=wmm_ca4.phcr, args=list(interval=rpyi, fmsy_proxy="f0.1", blim_proxy="f0.1", blim_proxy_multiplier = 0.25, nyears = 5)),
  # hcr: westmed MAP
  hcr = mseCtrl(method=westmedmap.hcr),
  # (i)mplementation (sys)tem: effort (in fact F if relationship E~F is 1)
  isys = mseCtrl(method=f.is)
)

nep9.rule <- mpCtrl(
  # estimation method: full feedback (a4a) sca
  est = mseCtrl(method=nep9.sa, args=list(fmodel=fmodel(nep9.fit),
                                          qmodel=lapply(qmodel(nep9.fit),"formula"), 
                                          srmodel=srmodel(nep9.fit), dfm=nep9.dfm, update=TRUE)),
  # parametrizing the HCR
  phcr = mseCtrl(method=wmm_ca4_bp.phcr, args=list(interval=rpyi, fmsy_proxy="f0.1")),
  # hcr: westmed MAP
  hcr = mseCtrl(method=westmedmap.hcr),
  # (i)mplementation (sys)tem: effort (in fact F if relationship E~F is 1)
  isys = mseCtrl(method=f.is)
)

nep11.rule <- mpCtrl(
  # estimation method: full feedback (a4a) sca
  est = mseCtrl(method=nep11.sa, args=list(fmodel=fmodel(nep11.fit),
                                           qmodel=lapply(qmodel(nep11.fit),"formula"), 
                                           srmodel=srmodel(nep11.fit), dfm=nep11.dfm, update=TRUE)),
  # parametrizing the HCR
  phcr = mseCtrl(method=wmm_ca4.phcr, args=list(interval=rpyi, fmsy_proxy="f0.1", blim_proxy="f0.1", blim_proxy_multiplier = 0.25)), #   # hcr: westmed MAP
  hcr = mseCtrl(method=westmedmap.hcr),
  # (i)mplementation (sys)tem: effort (in fact F if relationship E~F is 1)
  isys = mseCtrl(method=f.is)
)
#====================================================================
# Run simulations
#====================================================================
# note the way mp is coded it doesn't generate observations for the
# most recent year and it only generates observations for dy
loop_end <- fy-1-dl-ml
# restrictions for F
max_mult <- 2
min_mult <- 0.1

for(i in iy:(loop_end)){

  mseargs <- list(iy=i, fy=i+1, data_lag=dl, management_lag=ml, frq=af, real_iy=iy)

  if(exists("hke811.trck")) mseargs$mrun_track <- hke811.trck
  hke811.rule$isys@args$fmultiplier <- NULL
  hke811.mp <- mp(hke811.om, hke811.oem, hke811.iem, ctrl=hke811.rule, args=mseargs)

  if(exists("dps811.trck")) mseargs$mrun_track <- dps811.trck
  dps811.rule$isys@args$fmultiplier <- NULL
  dps811.mp <- mp(dps811.om, dps811.oem, dps811.iem, ctrl=dps811.rule, args=mseargs)

  if(exists("mut9.trck")) mseargs$mrun_track <- mut9.trck
  mut9.rule$isys@args$fmultiplier <- NULL
  mut9.mp <- mp(mut9.om, mut9.oem, mut9.iem, ctrl=mut9.rule, args=mseargs)
  
  if(exists("nep9.trck")) mseargs$mrun_track <- nep9.trck
  nep9.rule$isys@args$fmultiplier <- NULL
  nep9.mp <- mp(nep9.om, nep9.oem, nep9.iem, ctrl=nep9.rule, args=mseargs)
  
  if(exists("nep11.trck")) mseargs$mrun_track <- nep11.trck
  nep11.rule$isys@args$fmultiplier <- NULL
  nep11.mp <- mp(nep11.om, nep11.oem, nep11.iem, ctrl=nep11.rule, args=mseargs)

  #------------------------------------------------------------------
  # extract info to apply the multi stock HCR
  if(i==iy){
    hke811.trck0 <- tracking(hke811.mp)[year==ac(i)]
    dps811.trck0 <- tracking(dps811.mp)[year==ac(i)]
    mut9.trck0 <- tracking(mut9.mp)[year==ac(i)]
    nep9.trck0 <- tracking(nep9.mp)[year==ac(i)]
    nep11.trck0 <- tracking(nep11.mp)[year==ac(i)]
    } else {
    hke811.trck0 <- rbind(hke811.trck0, tracking(hke811.mp)[year==ac(i)])
    dps811.trck0 <- rbind(dps811.trck0, tracking(dps811.mp)[year==ac(i)])
    mut9.trck0 <- rbind(mut9.trck0, tracking(mut9.mp)[year==ac(i)])
    nep9.trck0 <- rbind(nep9.trck0, tracking(nep9.mp)[year==ac(i)])
    nep11.trck0 <- rbind(nep11.trck0, tracking(nep11.mp)[year==ac(i)])
  }

  #--------------------------------------------------------------------
  # multi HCR

  hke811.fstatus <- hke811.trck0[year %in% i & metric %in% c("F.est", "fmsy")]
  hke811.fstatus <- dcast(hke811.fstatus, year + iter ~ metric, value.var = "data")
  hke811.fstatus[, F_ratio := F.est / fmsy]

  mut9.fstatus <- mut9.trck0[year %in% i & metric %in% c("F.est", "fmsy")]
  mut9.fstatus <- dcast(mut9.fstatus, year + iter ~ metric, value.var = "data")
  mut9.fstatus[, F_ratio := F.est / fmsy]

  dps811.fstatus <- dps811.trck0[year %in% i & metric %in% c("F.est", "fmsy")]
  dps811.fstatus <- dcast(dps811.fstatus, year + iter ~ metric, value.var = "data")
  dps811.fstatus[, F_ratio := F.est / fmsy]
  
  nep9.fstatus <- nep9.trck0[year %in% i & metric %in% c("F.est", "fmsy")]
  nep9.fstatus <- dcast(nep9.fstatus, year + iter ~ metric, value.var = "data")
  nep9.fstatus[, F_ratio := F.est / fmsy]
  
  nep11.fstatus <- nep11.trck0[year %in% i & metric %in% c("F.est", "fmsy")]
  nep11.fstatus <- dcast(nep11.fstatus, year + iter ~ metric, value.var = "data")
  nep11.fstatus[, F_ratio := F.est / fmsy]

  # most vulnerable stock f multiplier
  fmult <- pmin(1/hke811.fstatus[,F_ratio], 1/dps811.fstatus[,F_ratio], 1/mut9.fstatus[,F_ratio], 1/nep9.fstatus[,F_ratio], 1/nep11.fstatus[,F_ratio])

  # average
#  fmult <- (1/hke811.fstatus[,F_ratio]+1/dps811.fstatus[,F_ratio]+1/mut9.fstatus[,F_ratio]+1/nep9.fstatus[,F_ratio]+1/nep11.fstatus[,F_ratio])/5
  
  # most representative stock f multiplier
#  fmult <- 1/hke811.fstatus[,F_ratio]

  # weighted average
  # hke811.cw <- hke811.trck0[year %in% i & metric == "C.est", "data"]
  # mut9.cw <- mut9.trck0[year %in% i & metric %in% c("C.est"), "data"]
  # dps811.cw <- dps811.trck0[year %in% i & metric %in% c("C.est"), "data"]
  # nep9.cw <- nep9.trck0[year %in% i & metric %in% c("C.est"), "data"]
  # nep11.cw <- nep11.trck0[year %in% i & metric %in% c("C.est"), "data"]
  # 
  # cw <- hke811.cw + mut9.cw + dps811.cw + nep9.cw + nep11.cw
  # 
  # hke811.cw <- hke811.cw/cw
  # mut9.cw <- mut9.cw/cw
  # dps811.cw <- dps811.cw/cw
  # nep9.cw <- nep9.cw/cw
  # nep11.cw <- nep11.cw/cw
  # 
  # fmult <- 1/hke811.fstatus[,F_ratio]*hke811.cw + 1/dps811.fstatus[,F_ratio]*dps811.cw +1/mut9.fstatus[,F_ratio]*mut9.cw+1/nep9.fstatus[,F_ratio]*nep9.cw+ 1/nep11.fstatus[,F_ratio]*nep11.cw

  #------------------------------------------------------------------
  # restricions for fmult
  fmult[fmult > max_mult] <- max_mult
  fmult[fmult < min_mult] <- min_mult
  
  ## prevent NaNs
  if(is.list(fmult)){
    fmult$data[!is.finite(fmult$data)] <- 1 # just for the first year
  } else{
    fmult[!is.finite(fmult)] <- 1
  }
  
  #------------------------------------------------------------------
  # project  
  hke811.rule$isys@args$fmultiplier <- unlist(fmult)
  dps811.rule$isys@args$fmultiplier <- unlist(fmult)
  mut9.rule$isys@args$fmultiplier <- unlist(fmult)
  nep9.rule$isys@args$fmultiplier <- unlist(fmult)
  nep11.rule$isys@args$fmultiplier <- unlist(fmult)
  mseargs <- list(iy=i, fy=i+dl+ml+1, data_lag=dl, management_lag=ml, frq=af, real_iy=iy)

  if(exists("hke811.trck")) mseargs$mrun_track <- hke811.trck
  hke811.mp <- mp(hke811.om, hke811.oem, hke811.iem, ctrl=hke811.rule, args=mseargs)
  hke811.om@stock[,ac(i+1)] <- hke811.mp@om@stock[,ac(i+1)]
  hke811.oem@observations[[1]][[1]][,ac(i+1)] <- hke811.mp@oem@observations[[1]][[1]][,ac(i+1)]
  hke811.oem@observations[[2]][,ac(i+1)] <- hke811.mp@oem@observations[[2]][,ac(i+1)]

  if(exists("dps811.trck")) mseargs$mrun_track <- dps811.trck
  dps811.mp <- mp(dps811.om, dps811.oem, dps811.iem, ctrl=dps811.rule, args=mseargs)
  dps811.om@stock[,ac(i+1)] <- dps811.mp@om@stock[,ac(i+1)]
  dps811.oem@observations[[1]][[1]][,ac(i+1)] <- dps811.mp@oem@observations[[1]][[1]][,ac(i+1)]
  dps811.oem@observations[[2]][,ac(i+1)] <- dps811.mp@oem@observations[[2]][,ac(i+1)]

  if(exists("mut9.trck")) mseargs$mrun_track <- mut9.trck
  mut9.mp <- mp(mut9.om, mut9.oem, mut9.iem, ctrl=mut9.rule, args=mseargs)
  mut9.om@stock[,ac(i+1)] <- mut9.mp@om@stock[,ac(i+1)]
  mut9.oem@observations[[1]][[1]][,ac(i+1)] <- mut9.mp@oem@observations[[1]][[1]][,ac(i+1)]
  mut9.oem@observations[[2]][,ac(i+1)] <- mut9.mp@oem@observations[[2]][,ac(i+1)]
  
  if(exists("nep9.trck")) mseargs$mrun_track <- nep9.trck
  nep9.mp <- mp(nep9.om, nep9.oem, nep9.iem, ctrl=nep9.rule, args=mseargs)
  nep9.om@stock[,ac(i+1)] <- nep9.mp@om@stock[,ac(i+1)]
  nep9.oem@observations[[1]][[1]][,ac(i+1)] <- nep9.mp@oem@observations[[1]][[1]][,ac(i+1)]
  nep9.oem@observations[[2]][,ac(i+1)] <- nep9.mp@oem@observations[[2]][,ac(i+1)]
  
  if(exists("nep11.trck")) mseargs$mrun_track <- nep11.trck
  nep11.mp <- mp(nep11.om, nep11.oem, nep11.iem, ctrl=nep11.rule, args=mseargs)
  nep11.om@stock[,ac(i+1)] <- nep11.mp@om@stock[,ac(i+1)]
  nep11.oem@observations[[1]][[1]][,ac(i+1)] <- nep11.mp@oem@observations[[1]][[1]][,ac(i+1)]
  nep11.oem@observations[[2]][,ac(i+1)] <- nep11.mp@oem@observations[[2]][,ac(i+1)]
  #------------------------------------------------------------------
  # keep track
  if(i==iy){
    hke811.trck <- tracking(hke811.mp)[year==ac(i)]
    mut9.trck <- tracking(mut9.mp)[year==ac(i)]
    dps811.trck <- tracking(dps811.mp)[year==ac(i)]
    nep9.trck <- tracking(nep9.mp)[year==ac(i)]
    nep11.trck <- tracking(nep11.mp)[year==ac(i)]
    } else {
    hke811.trck <- rbind(hke811.trck, tracking(hke811.mp)[year==ac(i)])
    mut9.trck <- rbind(mut9.trck, tracking(mut9.mp)[year==ac(i)])
    dps811.trck <- rbind(dps811.trck, tracking(dps811.mp)[year==ac(i)])
    nep9.trck <- rbind(nep9.trck, tracking(nep9.mp)[year==ac(i)])
    nep11.trck <- rbind(nep11.trck, tracking(nep11.mp)[year==ac(i)])
  }
  if(i==(loop_end)){
    tracking(hke811.mp) <- hke811.trck
    tracking(mut9.mp) <- mut9.trck
    tracking(dps811.mp) <- dps811.trck
    tracking(nep9.mp) <- nep9.trck
    tracking(nep11.mp) <- nep11.trck
  }
}

save(hke811.mp, hke811.om, mut9.mp, mut9.om, dps811.mp, dps811.om, 
     nep9.mp, nep9.om, nep11.mp, nep11.om, file="../results/basecase_EMU2.RData")

# stks <- window(FLStocks(hke811=hke811.om@stock, dps811=dps811.om@stock, mut9=mut9.om@stock), end=fy-1-dl-ml)
#
#  plot(stks) +
#    facet_wrap(qname~stock, ncol=3, scales = "free_y") +
#    scale_color_grey() +
#    scale_fill_grey() +
#    theme_bw()
#
# df0 <- data.frame(tracking(hke811.mp)[metric=="fmsy"])
# bwplot(data~year, data=df0)
#
# df0 <- data.frame(tracking(mut9.mp)[metric=="fmsy"])
# bwplot(data~year, data=df0)
#
# df0 <- data.frame(tracking(dps811.mp)[metric=="fmsy"])
# bwplot(data~year, data=df0)
#
# save.image(file="../results/MSEmulti250iter_15y_20262705.Rdata")
#
# hke811.fstatus <- tracking(hke811.mp)[metric %in% c("F.est", "fmsy")]
# hke811.fstatus <- dcast(hke811.fstatus, year + iter ~ metric, value.var = "data")
# hke811.fstatus[, F_ratio := F.est / fmsy]
# hke811.fstatus <- hke811.fstatus[, .(data = median(F_ratio, na.rm = TRUE)), by = year]
#
# mut9.fstatus <- tracking(mut9.mp)[metric %in% c("F.est", "fmsy")]
# mut9.fstatus <- dcast(mut9.fstatus, year + iter ~ metric, value.var = "data")
# mut9.fstatus[, F_ratio := F.est / fmsy]
# mut9.fstatus <- mut9.fstatus[, .(data = median(F_ratio, na.rm = TRUE)), by = year]
#
# dps811.fstatus <- tracking(dps811.mp)[metric %in% c("F.est", "fmsy")]
# dps811.fstatus <- dcast(dps811.fstatus, year + iter ~ metric, value.var = "data")
# dps811.fstatus[, F_ratio := F.est / fmsy]
# dps811.fstatus <- dps811.fstatus[, .(data = median(F_ratio, na.rm = TRUE)), by = year]
