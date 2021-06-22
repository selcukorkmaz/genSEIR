
#' Get COVID-19 Data
#'
#' The function collects the updated  COVID-19 data from the
#' John Hopkins University.

#' @param country name of the country. It should be a character string.
#' @param start a start date in mm/dd/yy format. Start date can not be earlier than 01/22/20. Start date can not be later than finish date. If start date is \code{NULL} then start date will be 01/22/20.
#' @param finish a finish date in mm/dd/yy format. Finish date can not be earlier than start date. If finish date is \code{NULL} then finish date will be the latest date at John-Hopkins CSSE system.
#'
#' @export getDataCOVID
#'
#' @importFrom utils read.csv
#'
#' @author Selcuk Korkmaz, \email{selcukorkmaz@gmail.com}
#'
#' @return a list of COVID-19 historical data including confirmed, death and recovered cases in desired time ranges.
#'
#'
#' @examples
#' covidData = getDataCOVID(country = "Italy",
#'                          start = "05/01/20",
#'                          finish = "12/31/20")
#' recovered = covidData$tableRecovered
#' deaths = covidData$tableDeaths
#' confirmed = covidData$tableConfirmed
#'
#' @seealso \code{\link{SEIQRDP}} \code{\link{fit_SEIQRDP}}
#'
#' @references Peng, L., Yang, W., Zhang, D., Zhuge, C., Hong, L. 2020. “Epidemic analysis of COVID-19 in China by dynamical modeling”, arXiv preprint arXiv:2002.06563.
#' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}


getDataCOVID <- function(country, start = NULL, finish = NULL){

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

  covidData = list(tableConfirmed = tableConfirmed, tableRecovered = tableRecovered, tableDeaths = tableDeaths)
  return(covidData)

}
