

tt<-SP_daily_asset_returns
tt<-tt[,2:ncol(tt)]


###########

ns <- duplicated(tt)
SP_daily_asset_returns <-
  SP_daily_asset_returns[which(!ns),]

SP_daily_asset_returns <-
  SP_daily_asset_returns[c(1,(3:nrow(SP_daily_asset_returns))),]


row.names(SP_daily_asset_returns) <-NULL


#########

save(SP_daily_asset_returns, file = 'SP_daily_asset_returns.RData', compress = 'xz')
