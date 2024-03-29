#' @importFrom stats frequency
SubATA.Damped <- function(train_set, pb, qb, model.Type, accuracy.Type, level.fix, trend.fix, trend.Search, phiStart, phiEnd, phiSize,
                          initialLevel, initialTrend, main_set, Holdout, HoldoutSet, Adjusted_P, h, Holdin, nmse, seas_periods, holdout_onestep)
{
  Xdata <- as.numeric(train_set)
  TA_0 <- Xdata-ATA.Shift(Xdata,1)
  TM_0 <- Xdata/ATA.Shift(Xdata,1)
  model.Type <- ifelse(is.null(model.Type),"B",model.Type)
  if (Holdout==TRUE){
    output <- SubATADampedHoldout(as.double(Xdata)
                                   , as.integer(ifelse(pb=="opt", -1, pb))
                                   , as.integer(ifelse(qb=="opt", -1, qb))
                                   , as.integer(switch(model.Type,"B"=0,"A"=1,"M"=2))
                                   , as.integer(switch(accuracy.Type,"MAE"=1,"MdAE"=2,"MSE"=3,"MdSE"=4,"MPE"=5,"MdPE"=6,"MAPE"=7,"MdAPE"=8,"sMAPE"=9,"sMdAPE"=10,"RMSE"=11,"MASE"=12,"OWA"=13,"AMSE"=14,"lik"=15,"sigma"=16,"GAMSE"=17))
                                   , as.integer(ifelse(level.fix, 1, 0))
                                   , as.integer(ifelse(trend.fix, 1, 0))
                                   , as.integer(ifelse(trend.Search, 1, 0))
                                   , as.double(phiStart)
                                   , as.double(phiEnd)
                                   , as.double(phiSize)
                                   , as.integer(switch(initialLevel,"none"=0,"mean"=1,"median"=2))
                                   , as.integer(switch(initialTrend,"none"=0,"mean"=1,"median"=2))
                                   , as.double(TA_0)
                                   , as.double(TM_0)
                                   , as.integer(seas_periods)
                                   , as.double(HoldoutSet)
                                   , as.integer(holdout_onestep))
  }else if (Holdin==TRUE){
    output <- SubATADampedHoldhin(as.double(Xdata)
                                   , as.integer(ifelse(pb=="opt", -1, pb))
                                   , as.integer(ifelse(qb=="opt", -1, qb))
                                   , as.integer(switch(model.Type,"B"=0,"A"=1,"M"=2))
                                   , as.integer(switch(accuracy.Type,"MAE"=1,"MdAE"=2,"MSE"=3,"MdSE"=4,"MPE"=5,"MdPE"=6,"MAPE"=7,"MdAPE"=8,"sMAPE"=9,"sMdAPE"=10,"RMSE"=11,"MASE"=12,"OWA"=13,"AMSE"=14,"lik"=15,"sigma"=16,"GAMSE"=17))
                                   , as.integer(ifelse(level.fix, 1, 0))
                                   , as.integer(ifelse(trend.fix, 1, 0))
                                   , as.integer(ifelse(trend.Search, 1, 0))
                                   , as.double(phiStart)
                                   , as.double(phiEnd)
                                   , as.double(phiSize)
                                   , as.integer(switch(initialLevel,"none"=0,"mean"=1,"median"=2))
                                   , as.integer(switch(initialTrend,"none"=0,"mean"=1,"median"=2))
                                   , as.double(TA_0)
                                   , as.double(TM_0)
                                   , as.integer(seas_periods)
                                   , as.integer(h)
                                   , as.integer(nmse))
  }else {
    output <- SubATADamped(as.double(Xdata)
                            , as.integer(ifelse(pb=="opt", -1, pb))
                            , as.integer(ifelse(qb=="opt", -1, qb))
                            , as.integer(switch(model.Type,"B"=0,"A"=1,"M"=2))
                            , as.integer(switch(accuracy.Type,"MAE"=1,"MdAE"=2,"MSE"=3,"MdSE"=4,"MPE"=5,"MdPE"=6,"MAPE"=7,"MdAPE"=8,"sMAPE"=9,"sMdAPE"=10,"RMSE"=11,"MASE"=12,"OWA"=13,"AMSE"=14,"lik"=15,"sigma"=16,"GAMSE"=17))
                            , as.integer(ifelse(level.fix, 1, 0))
                            , as.integer(ifelse(trend.fix, 1, 0))
                            , as.integer(ifelse(trend.Search, 1, 0))
                            , as.double(phiStart)
                            , as.double(phiEnd)
                            , as.double(phiSize)
                            , as.integer(switch(initialLevel,"none"=0,"mean"=1,"median"=2))
                            , as.integer(switch(initialTrend,"none"=0,"mean"=1,"median"=2))
                            , as.double(TA_0)
                            , as.double(TM_0)
                            , as.integer(seas_periods)
                            , as.integer(nmse))
  }
  ifelse(Holdout==TRUE & Adjusted_P==TRUE, new_pk <- round((output[1] * length(main_set))/ length(train_set)), new_pk <- output[1])
  ATA.last <- ATA.Core(main_set, pk = new_pk, qk = output[2], phik = output[3], mdlType = ifelse(output[4]==1,"A","M"), initialLevel = initialLevel, initialTrend = initialTrend)
  ATA.last$holdout <- Holdout
  ATA.last$holdin <- Holdin
  if(Holdout==TRUE){
    ATA.last$holdout.accuracy <- output[5]
    ATA.last$holdout.forecast <- ATAHoldoutForecast(as.double(Xdata)
                                                    , as.integer(output[1])
                                                    , as.integer(output[2])
                                                    , as.double(output[3])
                                                    , as.integer(output[4])
                                                    , as.integer(switch(initialLevel,"none"=0,"mean"=1,"median"=2))
                                                    , as.integer(switch(initialTrend,"none"=0,"mean"=1,"median"=2))
                                                    , as.double(TA_0)
                                                    , as.double(TM_0)
                                                    , as.integer(frequency(train_set))
                                                    , as.double(HoldoutSet)
                                                    , as.integer(holdout_onestep))
  }
  return(ATA.last)
}
