
# ------------------------------------------------------------
# Script:   0_postFilters.R
# Purpose:  Check similarity among valid names of marine species
# Updated:  15-06-2026
# ------------------------------------------------------------


#. Load packages ----

# Non-phonetic post-matching filters as implemented in Taxamatch
# https://doi.org/10.1371/journal.pone.0107510


## Step 3: Genus post-filter
genus_postFilter <- function(genus, preGenera, GeneraDist)
{
  # If distance exceed 3, discard the suggestion
  if(GeneraDist > 3)
  {
    return(F)
  }
  
  # Else...
  
  # Condition of 'good' characters (at least 50%)
  lpre <- min(c(nchar(genus), nchar(preGenera)))
  good <- (GeneraDist/lpre) <= 0.5
  
  # Condition of equal first character
  Tf1 <- substr(genus, 1, 1)
  If1 <- substr(preGenera, 1, 1)
  f1 <- Tf1==If1
  
  res <- F
  
  # If dist = 0 (exact match), keep match
  if(GeneraDist == 0)
  {
    res <- T
  }
  else if(good)
  {
    # If > 50% of characters match AND
    if(GeneraDist == 1)
    {
      # dist = 1, keep match
      res <- T
    }
    else if(GeneraDist < 4 & f1)
    {
      # Else, if dist = 2 OR 3 AND same first letter, keep match
      res <- T
    }
  }
  
  return(res)
}


## Step 6: Species post-filter
epithet_postFilter <- function(epithet, preEpithet, EpithetDist, GeneraDist)
{
  # If combined distance exceed 4, discard the suggestion
  if(EpithetDist + GeneraDist > 4 | EpithetDist > 4)
  {
    return(F)
  }
  
  # Else...
  
  # Condition of 'good' characters (at least 50%)
  lpre <- min(c(nchar(epithet), nchar(preEpithet)))
  good <- (EpithetDist/lpre) <= 0.5
  
  # Condition of equal firsts characters
  Tf1 <- substr(epithet, 1, 1)
  If1 <- substr(preEpithet, 1, 1)
  f1 <- Tf1==If1
  
  Tf3 <- substr(epithet, 1, 3)
  If3 <- substr(preEpithet, 1, 3)
  f3 <- Tf3==If3
  
  res <- F
  
  # If dist = 0 (exact match), keep match
  if(EpithetDist == 0)
  {
    res <- T
  }
  else if(good)
  {
    # If > 50% of characters match AND
    if(EpithetDist == 1)
    {
      # dist = 1, keep match
      res <- T
    }
    else if(EpithetDist %in% c(2,3) & f1)
    {
      # Else, if dist = 2 OR 3 AND same first letter, keep match
      res <- T
    }
    else if(EpithetDist == 4 & f3)
    {
      # Else, if dist = 4 AND same three first letters, keep match
      res <- T
    }
  }
  
  return(res)
}

