getDataCOVID <- function(start = NULL, finish = NULL, country = NULL){

  #' Get COVID-19 Data
  #'
  #' The function collects the updated  COVID-19 data from the
  #' John Hopkins University

  #' @param start enter a start date in mm/dd/yy format. Start date can not be earlier than 01/22/20. Start date can not be later than finish date. If start date is null then start date will be 01/22/20.
  #' @param finish enter a finish date in mm/dd/yy format. Finish date can not be earlier than start date. If finish date is null then finish date will be the latest date at John-Hopkins CSSE system.
  #' @param country name of the desired country.
  #'
  #' @export getDataCOVID
  #'
  #'
  #' @author Selçuk Korkmaz, \email{selcukorkmaz@gmail.com}
  #'
  #' @seealso \code{\link{SEIQRDP}} \code{\link{fit_SEIQRDP}}
  #'
  #' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}

  initialDate = as.Date("01/22/20", format = "%m/%d/%y")
  finalDate = as.Date(Sys.Date())-1
  time = seq(initialDate,finalDate, by="1 day")

  if(!is.null(start)){
    start = as.Date(start, format = "%m/%d/%y")
  }else{start = initialDate}

  if(!is.null(finish)){
    finish = as.Date(finish, format = "%m/%d/%y")
  }else{finish = finalDate}

  if(start > finish){
    stop(paste0("Start date is ", as.Date(start, format = "%m/%d/%y"), " and finish date is ", as.Date(finish, format = "%m/%d/%y"), ". Start date cannot be later than finish date."))
  }

  if(start < initialDate){
    stop(paste0("Start date is ", start, ". Start date cannot be earlier than ", initialDate))
  }

  if(finish > finalDate){
    stop(paste0("Finish date is ", finish, ". Finish date cannot be later than ", finalDate))
  }


  timeDF = cbind.data.frame(time,1:length(time))
  indT = timeDF[time >= start & time <= finish ,2]


  url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/"

  confirmed = read.csv(paste0(url, "time_series_covid19_confirmed_global.csv"), header = TRUE, check.names = FALSE)
  colnames(confirmed)[1:2] = c("ProvinceState", "CountryRegion")

  recovered = read.csv(paste0(url, "time_series_covid19_recovered_global.csv"), header = TRUE, check.names = FALSE)
  colnames(recovered)[1:2] = c("ProvinceState", "CountryRegion")

  deaths = read.csv(paste0(url, "time_series_covid19_deaths_global.csv"), header = TRUE, check.names = FALSE)
  colnames(deaths)[1:2] = c("ProvinceState", "CountryRegion")

  if(!is.null(country)){

    tableConfirmed = confirmed[confirmed$CountryRegion %in% country, c(1:4,indT+4)]
    tableRecovered = recovered[recovered$CountryRegion %in% country,c(1:4,indT+4)]
    tableDeaths = deaths[deaths$CountryRegion %in% country,c(1:4,indT+4)]

  }else{

    tableConfirmed = confirmed[, c(1:4,indT+4)]
    tableRecovered = recovered[,c(1:4,indT+4)]
    tableDeaths = deaths[,c(1:4,indT+4)]

  }

  return(list(tableConfirmed = tableConfirmed, tableRecovered=tableRecovered, tableDeaths=tableDeaths))

}
