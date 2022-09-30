#!/usr/bin/env ruby
#------------------------------------------------------------------------------
# parse log files and create an ascii file with the GRID job timing data 
# to be processed by su2020/scripts/grid_time_ana.C
#
# example: 
# --------
# grim/scripts/parse_grid_logs.rb -p ts_warm_bore -d 711_1010 -s s1 -j job [ --fileset=000] 
#
# output os stored in the xxx/timing_data subdirectory , at the same level as xxx/log
#
# comment: a bit kludgy, at this point, but works
#------------------------------------------------------------------------------
# puts "starting---"

require 'find'
require 'fileutils'
require 'getoptlong'

# puts " emoe"

$input_file = nil
$verbose    = nil


#-----------------------------------------------------------------------
def usage
  puts "usage: parse_grid_logs -d dsid [-v] "
  exit(-1)
end
#------------------------------------------------------------------------------
# specify defaults for the global variables and parse command line options
#------------------------------------------------------------------------------
def parse_command_line()
  opts = GetoptLong.new(
                        [ "--input-file"    , "-i",     GetoptLong::REQUIRED_ARGUMENT ],
                        [ "--verbose"       , "-v",     GetoptLong::NO_ARGUMENT       ]
                        )
  #----------------------------- defaults

  opts.each do |opt, arg|
    if    (opt == "--input-file"    ) ; $input_file = arg
    elsif (opt == "--verbose"       ) ; $verbose    = 1
    end

    if ($verbose != 0) ; puts "Option: #{opt}, arg #{arg.inspect}" ; end
  end
end


#------------------------------------------------------------------------------
def run(fn)

  puts "run : opening .#{fn}."
  f = File.open(fn);

  event_number = -1
  new_event = 0;
  e_gaas = 0;
  e_pd   = 0;


  puts "event/I:egaas/F:epd/F"

  f.each_line { |line|
    # if ($verbose) then puts "line=#{line}" end
    words = line.split();
    if (line.index("Begin processing")) then
      # new event
      if (new_event == 1) then
        if (e_gaas > 0) or (e_pd > 0) then
          puts "#{event_number} #{e_gaas} #{e_pd}"
        end
      end
      e_gaas       = 0;
      e_pd         = 0;
      new_event    = 1
      event_number = words[3].sub('th','').sub('nd','').sub('rd','').sub('st','');
    else 
      nw    = words.length
      if (nw > 9) then
        if    (words[8].index("GaasSensor")) then
          if ($verbose) then puts "line=#{line}" end
          e_gaas = e_gaas+words[5].to_f;
        elsif (words[8].index("InGaAs")) then
          if ($verbose) then puts "line=#{line}" end
          e_pd = e_pd+words[5].to_f;
        end
      end
    end
  }
  
  if (e_gaas > 0) or (e_pd > 0) then
    puts "#{event_number} #{e_gaas} #{e_pd}"
  end
  f.close
end

parse_command_line();

puts "input_file: #{$input_file}"


run($input_file);

exit(0)
