#'@title Plot credibility intervals from an \code{msocc} fit
#'
#'@description This function allows for visualization of credibility intervals
#'  for the derived parameters in an \code{msocc} fit. Optionally, one can check
#'  if each of the credibility intervals covers a value by specifying the
#'  \code{truth} option. This was originally designed for use in simulation,
#'  though could be used in other ways. \cr \cr By default, this funtion returns
#'  a list ggplots (the length of which depends upon whether \code{n} is
#'  specified).
#'
#'@param msocc an object of class \code{msocc}
#'@param level the level of the model to summarize; one of \code{c('site', 'sample', 'rep')}
#'@param truth optional vector of values to compare to each credibility interval; see details
#'@param n number of intervals to plot at once
#'@param quantiles quantiles to determine the level of credibility
#'@param burnin samples to discard as burn-in when summarizing the posterior
#'
#'@example examples/cred_plot_ex.R
#'
#'@return a \code{list} of ggplots
#'
#'@importFrom magrittr %>%
#'@export
#'
#'@md

cred_plot <- function(msocc, level = 'site', truth = NULL, n = 'all', quantiles = c(0.025, 0.975), burnin = 0){

  cred.lvl <- round(diff(quantiles) * 100)

  if(level == 'site'){
    if(n == 'all'){
      tmp <- posterior_summary(msocc, level = 'site', quantiles = quantiles, burnin = burnin)
      plot.df <- data.frame(
        site = tmp$site,
        mean = tmp$mean,
        lwr = tmp[,4],
        upr = tmp[,5]
      ) %>%
        dplyr::mutate(site = factor(site, levels = rev(site)))

      # ggthemr::ggthemr(palette = 'dust', layout = 'scientific')
      title <- paste0(cred.lvl, '% credibility intervals for psi by site')
      out <- plot.df %>%
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(x = lwr, y = site), shape = "|", size = 3, color = 'black') +
        ggplot2::geom_point(ggplot2::aes(x = upr, y = site), shape = "|", size = 3, color = 'black') +
        ggplot2::geom_segment(ggplot2::aes(x = lwr, y = site, xend = upr, yend = site), linetype = 'dashed') +
        ggplot2::geom_point(ggplot2::aes(y = site, x = mean), shape = 1) +
        ggplot2::labs(title = title,
             x = 'psi',
             y = 'site') +
        ggplot2::theme_bw()

      if(!is.null(truth)){
        truth.df <- data.frame(site = tmp$site, truth = truth)
        check.df <- data.frame(site = tmp$site, psi = -.05, col = ifelse(truth >= plot.df$lwr & truth <= plot.df$upr, 'green', 'red'))
        out <- out +
          ggplot2::geom_point(data = truth.df, ggplot2::aes(x = truth, y = site), shape = 16) +
          ggplot2::xlim(-0.05, 1) +
          ggplot2::geom_point(data = check.df, ggplot2::aes(x = psi, y = site), shape = 15, size = 2, color = check.df$col)
      } else{
        out <- out + ggplot2::xlim(0, 1)
      }
    } else{

      tmp <- posterior_summary(msocc, level = 'site', quantiles = quantiles, burnin = burnin)
      if(n > nrow(tmp)){
        plot.df <- data.frame(
          site = tmp$site,
          mean = tmp$mean,
          lwr = tmp[,4],
          upr = tmp[,5]
        )

        title <- paste0(cred.lvl, '% credibility intervals for psi by site')
        out <- plot.df %>%
          ggplot2::ggplot() +
          ggplot2::geom_point(ggplot2::aes(x = lwr, y = site), shape = "|", size = 3, color = 'black') +
          ggplot2::geom_point(ggplot2::aes(x = upr, y = site), shape = "|", size = 3, color = 'black') +
          ggplot2::geom_segment(ggplot2::aes(x = lwr, y = site, xend = upr, yend = site), linetype = 'dashed') +
          ggplot2::geom_point(ggplot2::aes(y = site, x = mean), shape = 1) +
          ggplot2::labs(title = title,
               x = 'psi',
               y = 'site') +
          ggplot2::theme_bw()

        if(!is.null(truth)){
          truth.df <- data.frame(site = tmp$site, truth = truth)
          check.df <- data.frame(site = tmp$site, psi = -.05, col = ifelse(truth >= plot.df$lwr & truth <= plot.df$upr, 'green', 'red'))
          out <- out +
            ggplot2::geom_point(data = truth.df, ggplot2::aes(x = truth, y = site), shape = 16) +
            ggplot2::xlim(-0.05, 1) +
            ggplot2::geom_point(data = check.df, ggplot2::aes(x = psi, y = site), shape = 15, size = 2, color = check.df$col)
        } else{
          out <- out + ggplot2::xlim(0, 1)
        }

      } else{

        if(!is.null(truth)){
          plot.df <- data.frame(
            site = tmp$site,
            mean = tmp$mean,
            lwr = tmp[,4],
            upr = tmp[,5],
            truth = truth
          ) %>%
            dplyr::mutate(site = factor(site, levels = rev(site)))
        } else{
          plot.df <- data.frame(
            site = tmp$site,
            mean = tmp$mean,
            lwr = tmp[,4],
            upr = tmp[,5]
          ) %>%
            dplyr::mutate(site = factor(site, levels = rev(site)))
        }

        split.df <- rep(1:ceiling(nrow(plot.df) / n), each = n)[1:nrow(plot.df)]
        df.list <- split(plot.df, split.df)

        title <- paste0(cred.lvl, '% credibility intervals for psi by site')
        out <- list()
        for(i in 1:length(df.list)){
          plot.df <- df.list[[i]]
          out[[i]] <- plot.df %>%
            ggplot2::ggplot() +
            ggplot2::geom_point(ggplot2::aes(x = lwr, y = site), shape = "|", size = 3, color = 'black') +
            ggplot2::geom_point(ggplot2::aes(x = upr, y = site), shape = "|", size = 3, color = 'black') +
            ggplot2::geom_segment(ggplot2::aes(x = lwr, y = site, xend = upr, yend = site), linetype = 'dashed') +
            ggplot2::geom_point(ggplot2::aes(y = site, x = mean), shape = 1) +
            ggplot2::labs(title = title,
                 x = 'psi',
                 y = 'site') +
            ggplot2::theme_bw()

          if(!is.null(truth)){
            check.df <- data.frame(site = plot.df$site, psi = -.05, col = ifelse(plot.df$truth >= plot.df$lwr & plot.df$truth <= plot.df$upr, 'green', 'red'))
            out[[i]] <- out[[i]] +
              ggplot2::geom_point(data = plot.df, ggplot2::aes(x = truth, y = site), shape = 16) +
              ggplot2::xlim(-0.05, 1) +
              ggplot2::geom_point(data = check.df, ggplot2::aes(x = psi, y = site), shape = 15, size = 2, color = check.df$col)
          } else{
            out[[i]] <- out[[i]] + ggplot2::xlim(0, 1)
          }
        }
      }
    }
  }

  if(level == 'sample'){

    if(n == 'all'){
      tmp <- posterior_summary(msocc, level = 'sample', quantiles = quantiles, burnin = burnin)
      plot.df <- data.frame(
        site_sample = paste(tmp$site, tmp$sample, sep = "_"),
        mean = tmp$mean,
        lwr = tmp[,6],
        upr = tmp[,7]
      ) %>%
        dplyr::mutate(site_sample = factor(site_sample, levels = rev(site_sample)))

      # ggthemr::ggthemr(palette = 'dust', layout = 'scientific')
      title <- paste0(cred.lvl, '% credibility intervals for theta by site-sample combination')
      out <- plot.df %>%
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(x = lwr, y = site_sample), shape = "|", size = 3, color = 'black') +
        ggplot2::geom_point(ggplot2::aes(x = upr, y = site_sample), shape = "|", size = 3, color = 'black') +
        ggplot2::geom_segment(ggplot2::aes(x = lwr, y = site_sample, xend = upr, yend = site_sample), linetype = 'dashed') +
        ggplot2::geom_point(ggplot2::aes(y = site_sample, x = mean), shape = 1) +
        ggplot2::labs(title = title,
             x = 'theta',
             y = 'site_sample') +
        ggplot2::theme_bw()

      if(!is.null(truth)){
        truth.df <- data.frame(site = plot.df$site_sample, truth = truth)
        check.df <- data.frame(site = plot.df$site_sample, theta = -.05, col = ifelse(truth >= plot.df$lwr & truth <= plot.df$upr, 'green', 'red'))
        out <- out +
          ggplot2::geom_point(data = truth.df, ggplot2::aes(x = truth, y = site), shape = 16) +
          ggplot2::xlim(-0.05, 1) +
          ggplot2::geom_point(data = check.df, ggplot2::aes(x = theta, y = site), shape = 15, size = 2, color = check.df$col)
      } else{
        out <- out + ggplot2::xlim(0,1)
      }
    } else{

      tmp <- posterior_summary(msocc, level = 'sample', quantiles = quantiles, burnin = burnin)
      if(n > nrow(tmp)){
        plot.df <- data.frame(
          site_sample = paste(tmp$site, tmp$sample, sep = "_"),
          mean = tmp$mean,
          lwr = tmp[,6],
          upr = tmp[,7]
        ) %>%
          dplyr::mutate(site_sample = factor(site_sample, levels = rev(site_sample)))

        # ggthemr::ggthemr(palette = 'dust', layout = 'scientific')
        title <- paste0(cred.lvl, '% credibility intervals for theta by site-sample combination')
        out <- plot.df %>%
          ggplot2::ggplot() +
          ggplot2::geom_point(ggplot2::aes(x = lwr, y = site_sample), shape = "|", size = 3, color = 'black') +
          ggplot2::geom_point(ggplot2::aes(x = upr, y = site_sample), shape = "|", size = 3, color = 'black') +
          ggplot2::geom_segment(ggplot2::aes(x = lwr, y = site_sample, xend = upr, yend = site_sample), linetype = 'dashed') +
          ggplot2::geom_point(ggplot2::aes(y = site_sample, x = mean), shape = 1) +
          ggplot2::labs(title = title,
               x = 'theta',
               y = 'site_sample') +
          ggplot2::theme_bw()

        if(!is.null(truth)){
          truth.df <- data.frame(site = plot.df$site_sample, truth = truth)
          check.df <- data.frame(site = plot.df$site_sample, theta = -.05, col = ifelse(truth >= plot.df$lwr & truth <= plot.df$upr, 'green', 'red'))
          out <- out +
            ggplot2::geom_point(data = truth.df, ggplot2::aes(x = truth, y = site), shape = 16) +
            ggplot2::xlim(-0.05, 1) +
            ggplot2::geom_point(data = check.df, ggplot2::aes(x = theta, y = site), shape = 15, size = 2, color = check.df$col)
        } else{
          out <- out + ggplot2::xlim(0,1)
        }
      } else{
        if(!is.null(truth)){
          plot.df <- data.frame(
            site_sample = paste(tmp$site, tmp$sample, sep = "_"),
            mean = tmp$mean,
            lwr = tmp[,6],
            upr = tmp[,7],
            truth = truth
          ) %>%
            dplyr::mutate(site_sample = factor(site_sample, levels = rev(site_sample)))
        } else{
          plot.df <- data.frame(
            site_sample = paste(tmp$site, tmp$sample, sep = "_"),
            mean = tmp$mean,
            lwr = tmp[,6],
            upr = tmp[,7]
          ) %>%
            dplyr::mutate(site_sample = factor(site_sample, levels = rev(site_sample)))
        }

        split.df <- rep(1:ceiling(nrow(plot.df) / n), each = n)[1:nrow(plot.df)]
        df.list <- split(plot.df, split.df)

        title <- paste0(cred.lvl, '% credibility intervals for theta by site-sample combination')
        out <- list()
        for(i in 1:length(df.list)){
          plot.df <- df.list[[i]]
          out[[i]] <- plot.df %>%
            ggplot2::ggplot() +
            ggplot2::geom_point(ggplot2::aes(x = lwr, y = site_sample), shape = "|", size = 3, color = 'black') +
            ggplot2::geom_point(ggplot2::aes(x = upr, y = site_sample), shape = "|", size = 3, color = 'black') +
            ggplot2::geom_segment(ggplot2::aes(x = lwr, y = site_sample, xend = upr, yend = site_sample), linetype = 'dashed') +
            ggplot2::geom_point(ggplot2::aes(y = site_sample, x = mean), shape = 1) +
            ggplot2::labs(title = title,
                 x = 'theta',
                 y = 'site_sample') +
            ggplot2::theme_bw()

          if(!is.null(truth)){
            check.df <- data.frame(site_sample = plot.df$site_sample, theta = -.05, col = ifelse(plot.df$truth >= plot.df$lwr & plot.df$truth <= plot.df$upr, 'green', 'red'))
            out[[i]] <- out[[i]] +
              ggplot2::geom_point(data = plot.df, ggplot2::aes(x = truth, y = site_sample), shape = 16) +
              ggplot2::xlim(-0.05, 1) +
              ggplot2::geom_point(data = check.df, ggplot2::aes(x = theta, y = site_sample), shape = 15, size = 2, color = check.df$col)
          } else{
            out[[i]] <- out[[i]] + ggplot2::xlim(0, 1)
          }
        }
      }
    }
  }

  if(level == 'rep'){

    if(n == 'all'){
      tmp <- posterior_summary(msocc, level = 'rep', quantiles = quantiles, burnin = burnin)
      plot.df <- data.frame(
        site_sample_rep = paste(tmp$site, tmp$sample, tmp$rep, sep = "_"),
        mean = tmp$mean,
        lwr = tmp[,6],
        upr = tmp[,7]
      ) %>%
        dplyr::mutate(site_sample_rep = factor(site_sample_rep, levels = rev(site_sample_rep)))

      # ggthemr::ggthemr(palette = 'dust', layout = 'scientific')
      title <- paste0(cred.lvl, '% credibility intervals for p by site-sample-rep combination')
      out <- plot.df %>%
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(x = lwr, y = site_sample_rep), shape = "|", size = 3, color = 'black') +
        ggplot2::geom_point(ggplot2::aes(x = upr, y = site_sample_rep), shape = "|", size = 3, color = 'black') +
        ggplot2::geom_segment(ggplot2::aes(x = lwr, y = site_sample_rep, xend = upr, yend = site_sample_rep), linetype = 'dashed') +
        ggplot2::geom_point(ggplot2::aes(y = site_sample_rep, x = mean), shape = 1) +
        ggplot2::labs(title = title,
             x = 'p',
             y = 'site_sample_rep') +
        ggplot2::theme_bw()

      if(!is.null(truth)){
        truth.df <- data.frame(site = plot.df$site_sample_rep, truth = truth)
        check.df <- data.frame(site = plot.df$site_sample_rep, p = -.05, col = ifelse(truth >= plot.df$lwr & truth <= plot.df$upr, 'green', 'red'))
        out <- out +
          ggplot2::geom_point(data = truth.df, ggplot2::aes(x = truth, y = site), shape = 16) +
          ggplot2::xlim(-0.05, 1) +
          ggplot2::geom_point(data = check.df, ggplot2::aes(x = p, y = site), shape = 15, size = 2, color = check.df$col)
      } else{
        out <- out + ggplot2::xlim(0,1)
      }
    } else{
      tmp <- posterior_summary(msocc, level = 'rep', quantiles = quantiles, burnin = burnin)
      if(n > nrow(tmp)){
        plot.df <- data.frame(
          site_sample_rep = paste(tmp$site, tmp$sample, tmp$rep, sep = "_"),
          mean = tmp$mean,
          lwr = tmp[,6],
          upr = tmp[,7]
        ) %>%
          dplyr::mutate(site_sample_rep = factor(site_sample_rep, levels = rev(site_sample_rep)))

        # ggthemr::ggthemr(palette = 'dust', layout = 'scientific')
        title <- paste0(cred.lvl, '% credibility intervals for p by site-sample-rep combination')
        out <- plot.df %>%
          ggplot2::ggplot() +
          ggplot2::geom_point(ggplot2::aes(x = lwr, y = site_sample_rep), shape = "|", size = 3, color = 'black') +
          ggplot2::geom_point(ggplot2::aes(x = upr, y = site_sample_rep), shape = "|", size = 3, color = 'black') +
          ggplot2::geom_segment(ggplot2::aes(x = lwr, y = site_sample_rep, xend = upr, yend = site_sample_rep), linetype = 'dashed') +
          ggplot2::geom_point(ggplot2::aes(y = site_sample_rep, x = mean), shape = 1) +
          ggplot2::labs(title = title,
               x = 'p',
               y = 'site_sample_rep') +
          ggplot2::theme_bw()

        if(!is.null(truth)){
          truth.df <- data.frame(site = plot.df$site_sample_rep, truth = truth)
          check.df <- data.frame(site = plot.df$site_sample_rep, p = -.05, col = ifelse(truth >= plot.df$lwr & truth <= plot.df$upr, 'green', 'red'))
          out <- out +
            ggplot2::geom_point(data = truth.df, ggplot2::aes(x = truth, y = site), shape = 16) +
            ggplot2::xlim(-0.05, 1) +
            ggplot2::geom_point(data = check.df, ggplot2::aes(x = p, y = site), shape = 15, size = 2, color = check.df$col)
        } else{
          out <- out + ggplot2::xlim(0,1)
        }
      } else{
        if(!is.null(truth)){
          plot.df <- data.frame(
            site_sample_rep = paste(tmp$site, tmp$sample, tmp$rep, sep = "_"),
            mean = tmp$mean,
            lwr = tmp[,6],
            upr = tmp[,7],
            truth = truth
          ) %>%
            dplyr::mutate(site_sample_rep = factor(site_sample_rep, levels = rev(site_sample_rep)))
        } else{
          plot.df <- data.frame(
            site_sample_rep = paste(tmp$site, tmp$sample, tmp$rep, sep = "_"),
            mean = tmp$mean,
            lwr = tmp[,6],
            upr = tmp[,7]
          ) %>%
            dplyr::mutate(site_sample_rep = factor(site_sample_rep, levels = rev(site_sample_rep)))
        }

        split.df <- rep(1:ceiling(nrow(plot.df) / n), each = n)[1:nrow(plot.df)]
        df.list <- split(plot.df, split.df)

        title <- paste0(cred.lvl, '% credibility intervals for p by site-sample-rep combination')
        out <- list()
        for(i in 1:length(df.list)){
          plot.df <- df.list[[i]]
          out[[i]] <- plot.df %>%
            ggplot2::ggplot() +
            ggplot2::geom_point(ggplot2::aes(x = lwr, y = site_sample_rep), shape = "|", size = 3, color = 'black') +
            ggplot2::geom_point(ggplot2::aes(x = upr, y = site_sample_rep), shape = "|", size = 3, color = 'black') +
            ggplot2::geom_segment(ggplot2::aes(x = lwr, y = site_sample_rep, xend = upr, yend = site_sample_rep), linetype = 'dashed') +
            ggplot2::geom_point(ggplot2::aes(y = site_sample_rep, x = mean), shape = 1) +
            ggplot2::labs(title = title,
                 x = 'p',
                 y = 'site_sample_rep') +
            ggplot2::theme_bw()

          if(!is.null(truth)){
            check.df <- data.frame(site_sample_rep = plot.df$site_sample_rep, p = -.05, col = ifelse(plot.df$truth >= plot.df$lwr & plot.df$truth <= plot.df$upr, 'green', 'red'))
            out[[i]] <- out[[i]] +
              ggplot2::geom_point(data = plot.df, ggplot2::aes(x = truth, y = site_sample_rep), shape = 16) +
              ggplot2::xlim(-0.05, 1) +
              ggplot2::geom_point(data = check.df, ggplot2::aes(x = p, y = site_sample_rep), shape = 15, size = 2, color = check.df$col)
          } else{
            out[[i]] <- out[[i]] + ggplot2::xlim(0, 1)
          }
        }
      }
    }
  }

  return(out)
}


