% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_vac_schedule.R
\name{convert_vac_schedule}
\alias{convert_vac_schedule}
\title{convert cumulative vaccination schedule to non-cumulative ----------------------------------}
\usage{
convert_vac_schedule(
  vac_schedule,
  ve,
  hosp_multiplier,
  delay,
  ve_trans,
  add_child_vac = FALSE,
  child_vac_coverage = 0.75,
  child_doses_per_day = 50000,
  child_vac_start_date = "2021-09-01",
  wane = FALSE,
  k = 0.03,
  t0 = 180,
  add_extra_dates = FALSE,
  extra_end_date = "2022-03-31"
)
}
\arguments{
\item{vac_schedule}{a data frame that the proportion of the population who receives vaccines
at each time point. Rows are time points, columns are vaccine type, age group, and dose. For example,
the first column is the proportion of individuals in age group 1 who receive dose 1 of the first vaccine
type. The function assumes 9 age groups and four vaccine types, each with a 2-dose regimen.}

\item{ve}{a named list of vaccine effectiveness against infection for each dose of each vaccine type.}

\item{hosp_multiplier}{a named list of the values for the vaccine effectiveness against hospitalization
for each dose of each vaccine type. Vaccine effectiveness against hospitalization is incorporated as a 
multiplier on the probability of being hospitalized after infection as (1 – VE_hospitalization) divided by (1-VE_infection)
to account for the the inclusion of people who are never infected (and thus never hospitalized) included 
in the estimation of VE against hospitalization.}

\item{delay}{a named list of the time to protection for each dose of each vaccine type.}

\item{ve_trans}{a named list of vaccine effectiveness against transmission for each dose of each vaccine type.}

\item{add_child_vac}{logical, if TRUE 5-11 year olds are vaccinated}

\item{child_vac_coverage}{total vaccine coverage in 5-11 to achieve (between 0 and 1) if add_child_vac = TRUE}

\item{child_doses_per_day}{the number of doses to be administered to 5-11 year olds if add_child_vac = TRUE}

\item{child_vac_start_date}{character string of the date (YYYY-MM-DD format) to start vaccinating 5-11 year
olds, if add_child_vac = TRUE}

\item{wane}{logical, if TRUE vaccine effectiveness wanes by a logistic function parameterized by arguments
k and t0.}

\item{k}{logistic growth rate}

\item{t0}{the time point at the midpoint of the logistic curve (where 50\% waning occurs)}

\item{add_extra_dates}{logical, if TRUE add extra rows to the vaccination schedule until extra_end_date}

\item{extra_end_date}{character string of the date (YYYY-MM-DD format) to end the vaccination schedule (which will
also be the last date of the simulation)}
}
\value{
list of vaccination rate by day and dose and weighted VE and delay to protection by day and dose
}
\description{
convert cumulative vaccination schedule to non-cumulative ----------------------------------
}
\keyword{vacamole}
