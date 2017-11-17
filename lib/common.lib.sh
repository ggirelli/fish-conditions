#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 2.0.0
# Date: 20171118
# Description: functions for common operations
# 
# ------------------------------------------------------------------------------



# FUNCTIONS ====================================================================

function confirm_settings() {
    opt_string=${1}

    settings_confirm=`echo -e "$opt_string" | sed -E 's/^/ /'`
    settings_confirm="
    ##############################################
    #                                            #
    #  PLEASE, DOUBLE CHECK THE SETTINGS BELOW   #
    # (press 'q' to continue when you are ready) #
    #                                            #
    ##############################################

    $settings_confirm"

    echo -e "$settings_confirm" | less

    msg="$msg\nRun the analysis?\nYes (y), Abort (a), Show again (s)"
    clear; echo -e $msg; read -e ans

    end=0
    while [[ 0 -eq $end ]]; do
        if [[ -z $ans ]]; then
          echo -e $msg
          read -e ans
        elif [[ 'a' == $ans ]]; then
          end=1
          echo "Aborted."
          exit 1
        elif [[ 'y' == $ans ]]; then
          echo -e "\n"
          end=1
        elif [[ 's' == $ans ]]; then
          echo -e "$settings_confirm" | less
          clear; echo -e $msg; read -e ans
        else
          echo -e $msg
          read -e ans
        fi
    done
}

function confirm_overwrite() {
    outdir=${1}

    if [ -e "$outdir" ]; then
        end=0
    else
        end=1
    fi
    while [[ 0 -eq end ]]; do
        if [ -d "$outdir" ]; then
            msg="The specified folder already exists."
            msg="$msg\nFolder: $outdir"
            msg="$msg\nOverwrite?\nYes (y), No (n)"
            echo -e $msg; read -e ans

            while [[ 0 -eq $end ]]; do
                if [[ -z $ans ]]; then
                  echo -e $msg
                  read -e ans
                elif [[ 'n' == $ans ]]; then
                  end=1
                  msg="Aborted. Please, provide a different output folder path."
                  echo $msg
                  exit 1
                elif [[ 'y' == $ans ]]; then
                  echo -e "\n"
                  rm -r "$outdir"
                  end=1
                else
                  echo -e $msg
                  read -e ans
                fi
            done
        else
            msg="The specified folder already exists as a file."
            msg="$msg\nFolder: $outdir"
            msg="$msg\nOverwrite?\nYes (y), No (n)"
            echo -e $msg; read -e ans

            end=0
            while [[ 0 -eq $end ]]; do
                if [[ -z $ans ]]; then
                  echo -e $msg
                  read -e ans
                elif [[ 'n' == $ans ]]; then
                  end=1
                  msg="Aborted. Please, provide a different output folder path."
                  echo "$msg"
                  exit 1
                elif [[ 'y' == $ans ]]; then
                  echo -e "\n"
                  rm "$outdir"
                  end=1
                else
                  echo -e $msg
                  read -e ans
                fi
            done
        fi
    done
}

################################################################################