function []= x_cluster()
% x_cluster.m    AS 2.2i                                                 10:29 18.11.2011
%                Cluster-Analysis including PCA, MDS and k-means-clustering
%-- { ------------------------------------------------------------------------------------
% INPUT-format of data (feature vectors, distance matrix, grand-mean reference data):
%    * Rows:    Cases. The dissimilarities between cases will be analysed (clustered)
%                      Not individual cases, but GROUP-MEANS
%      Columns: Case-names (col-1) or
%               parameters (col-2:end), e.g. receptor densities
%    * ASCII-format, according to the MatLab tblread()format:
%      <[data,varnames,casenames] = tblread(filename)>
%    * Accepts files in the current wd <.\>, or as indicated by the complete pathname.
%    * No grouping variable, 
%                          ** only group means are accepted **.
%      Means must be generated using an external program, eg in SYSTAT 
%      <DATA_Prepare_multivar_group-means_(syz-inp).syc>. No normalization, no z-sco. 
%    * Avoid leading spaces in the data of thre first row! See examples below
%    * The first column is always interpreted as case_names, even if numeric!
%      If no meaningfull names are available, use case numbers 1,2,3 ... or dummies
%    * Numeric data starts always at cell (2,2)
%
%    tblread(): Return-Value  Description
%    -----------------------------------------------------------------------------------
%    data          Numeric matrix with a value for each variable-case pair (cell).
%    varnames      String matrix containing the variable names,     the first row
%    casenames     String matrix containing the names of each case, the first column
%    -----------------------------------------------------------------------------------
%    Examples of input files:
%
%    Example 'feature_vector'. String 'SPECIES' will be ignored!:     -------------->
%    SPECIES     MEANX_O  SD_O    ...
%    alo_a17_1   0.1148  -0.1131  ...
%    aot_a17_1  -0.2947   0.0146  ...
%    ate_a17_1   1.0774  -0.2705  ...
%    bon_a17_1  -1.4928   0.8973  ...
%    cal_a17_1  -0.4304   0.6356  ...
%
%    Example 'distance matrix'. Complete header line will be ignored!:  ------------>
%    NAME      Area_1          Area_2          Area_3         ...
%    Area_1    0.0000000e+000  2.7256491e+000  1.7926109e+000 ...
%    Area_2    2.7256491e+000  0.0000000e+000  1.8544197e+000 ...
%    Area_3    1.7926109e+000  1.8544197e+000  0.0000000e+000 ...
%
%    Example 'grand-mean reference data': ------------------------------------------>
%      Demonstration of data used for normalization (optional):
%      String 'STATISTIC' is ignored, but at least one blank must start the line!
%    * Path to this file is fixed to [MLMNETPATH '\data\']
%    * Data gives means (averages, first line) and medians (secons line)
%
%      STATISTIC   AMPA   KAIN   MK80   MUSC   CGP5   FLUM ...
%      MEAN       558.3  570.5 1175.1 1555.1 2170.0 2552.4 ...
%      MEDIAN     533.0  521.0 1171.2 1520.0 2149.0 2581.0 ...

%----------------------------------------------------------------------------------------
% OUTPUT:
%    All files (data, plots) created via <x_cluster.m> will be stored to <.\OUT>.
%
%    'name': name of input-file read via tblread()                    ['a17_samples.txt']
%    * Input data in *.mat format:
%                                                 ['a17_sample.txt' -> 'a17_samples.mat']
%    * Normalized data, tbl-format. Example z-scores             ['a17_samples_zsco.csv']
%    * Graph PCA  (optional)                                        [a17_samples_pca.emf]
%    * Graph MDS  (optional)                                        [a17_samples_mds.emf]
%    * Dendrogram (optional)                                      [a17_samples_hclus.emf]
%    * Silhouett plot of k-means clustering  (optional)           [a17_samples_kclus.emf]
%    * Matrix of Euclidean distances tbl-format           [a17_samples_'metric'-dist.csv]
%      'metric' is euc/seu/cos/cor/cit (first 3 chars from <dist_metric>)
%                                                              [a17_samples_euc-dist.csv]
%----------------------------------------------------------------------------------------
% SUBFUNCTIONS
%     y_get_axes_limits()   Get axes limits from current figure
%     y_showgraph():        Interrupts program, shows graph and ask for action
%     y_tblread_block():    Reads input data in tblread() format.
%     z_setgraph():         Definition of Plotsize
%----------------------------------------------------------------------------------------
% ACTIONS:
%   * Normalize data (Opt.): none/divide-by-SD/z-scores/divide-by-MEAN
%   * Transpose data (Opt.): Not tested! For clustering of variables (feature_names)
%                            A clustering of cases (rows) is performed if  transpose is
%                            not activated
%   * Boxplot        (Opt.): Statistics of each column is plotted as a boxplot to detect
%----------------------------------------------------------------------------------------
% TEST-DATA:
%   area_means_only_no42.txt            Feature vectors: autoradiography   (LG)
%   digrju_dist_mat_euc1_cluster.txt    Distance matrix, 6 objects (areas) (JU)
%                                       tbl-format
%----------------------------------------------------------------------------------------
% BOXPLOT-interpretation
%  * Boxplot(X) produces a box and whisker plot for each column of X.
%    The box has lines at the lower quartile, median, and upper quartile values.
%    The whiskers are lines extending from each end of the box to show the extent of the
%    rest of the data. Outliers are data with values beyond the ends of the whiskers.
%    If there is no data outside the whisker, a dot is placed at the bottom whisker.
%----------------------------------------------------------------------------------------
% DATA used (tables, LUT's etc):
%  * Grand means in [MLMNETPATH '\data\']
%    grand_mean_cortex.csv  Grand mean cortex data, means and median
%----------------------------------------------------------------------------------------
% HISTORY:
%   04-12-13    MDS: non-linear/linear mds (ML7.1) including Shepard plot
%   05-10-05    Input data:
%               distance matrix, for MDS + Hierarchical-Cluster analysis only.
%   .....
%   06-09-20    Fixed a bug in the grand-mean/median normalization procedure
%   07-05-10    Output of matrix Euclidean distances (via hier-clust)
%   11-11-17    fileparts() adjusted to <7.12.0.635 (R2011a)>
%   12-08-08    corrected user entries for norm_mode:  0= no normalization, 1= divide by sd
%----------------------------------------------------------------------------------------
% COMMENTS:
%   * Linkage:
%     When linkage is 'centroid', 'median', or 'ward', the output of linkage is
%     meaningful only if the input Y contains Euclidean distances
%----------------------------------------------------------------------------------------
% BUGS AND LIMITATIONS:
%   *
%- } -------------------------------------------------------------------------------------

% { ------------------------------------------------------------------------------- m+main
global DIARY_NAME      % Name output file (via diary)
global GRAPH_FORMAT    % Fileformat for graphs (output)                             [emf]
global FONT_SIZE       % Fontsize for labelling graphes
global MLMNETPATH      % Path to external y_*.m and used data               [D:\MLMFILES]
                       % data (means and mediand for grand mean) in [MLMNETPATH '\data\']
global N_COLUMNS
global N_ROWS
global PATH_OUT        % path to output data            [E:\DTN\MATLAB_DEV\x_cluster\OUT]
global WORK_DIR        % Working directory                  [E:\DTN\MATLAB_DEV\x_cluster]

%-------------------------------------------------------------------------  m+globals-set
DIARY_NAME   = 'x_cluster.mof';
FONT_SIZE    = 8;   % to be set in m~user
GRAPH_FORMAT = 'emf';  % emf/tif/png
%--------------------------------------------------------------------------- m+initialize
clc          % Clear command window
close all    % Clear all figures

if (exist(DIARY_NAME,'file')); delete(DIARY_NAME); end
s_line=repmat('-',1,80);
diary(DIARY_NAME); diary on
% Diary header (program identification)
  s_text = [s_line '\nProgram:\tx_cluster.m' ...
                        '\nDate:   \t' datestr(now) '\n'];
  fprintf(1,s_text); fprintf(1,[s_line '\n']);

%-------------------------------------------------------------------------------- m+mkdir
% Create directory for output, fixed to [WORK_DIR '\OUT'], name is <PATH_OUT>
% WORK_DIR is the current working directory (.\)
  WORK_DIR = pwd;
  PATH_OUT= [WORK_DIR '\OUT']; flag_true = exist(PATH_OUT,'dir');
  if     (flag_true ); fprintf(1,'.... Directory %s found\n',PATH_OUT)
  else    mkdir(WORK_DIR,'\OUT'); fprintf(1,'.... Directory %s created\n',PATH_OUT)
  end

%>{>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> m+user #+mm
%  do_boxplot         'yes'/'no' Create boxplot?
%  do_scatter_plot    'no'/'pca'/'mds' Create 2D-Scatterplot?
%                     'pca'  Generate Principal Component scores (pca)
%                     'mds'  Multidimensional Scaling (mds), linear or nonlinear (default)
%  do_h_cluster       'yes'/'no' Apply hierarchical clustering with dendrogram?
%  do_kmeans_cluster  'yes'/'no' Perform k-means clustering?
%  data_mode          Type of input data:
%                     'feature_vector' / 'distance_matrix'
%  data_file          Name of ASCII-input file with extension (*.{txt,csv})
%                     Format must conform the tblread() format definition (see above)
%  col_varfirst       First column of the block of variables to be processed          [2]
%  col_varlast        Last  column of the block of variables to be processed          [8]
%  norm_mode          0/1/2 Normalization of input data (group means are normalized!)
%                     0/none; 1/divide by SD;  2/z-scores; 3/divide by sample mean
%                             4/divide by grand mean ;     5/divide by grand median
%  save_cluster       'yes'/'no' Save Cluster number (index) of each case ?        ['no']
%  dist_metric   >>   Distance metric used in H-Clust/mds/K-Means-Clust:
%                       'euclid'/'seuclid'/'cosine'/'correlation'/'cityblock'
%  link_method        'single'/'complete'/'average'/'centroid'/'ward'
%  mds_criterion      Multidimensional scaling:
%                     'stress'/'metricstress' --> for non-metric/metric scaling
%  maxn_numbers       In PCA, only fewer data points will be labeled                 [36]
%  FONT_SIZE          Font size to label data points and tick labels                  [8]
%------------------------------------------------------------------------- user entries >
% * Data input in tbl-format: 1/feature file, 2/grand mean reference data (optional):
%   Columns are counted with respect to the original data file including:
%   1st row:    variable name,
%   1st column: case names.     Numerical data starts at cell (2,2)
% --------------- Actions (only the first 2 chars are used)
do_boxplot        = 'no';               % 'yes'/'no'
do_scatter_plot   = 'pca';               % 'no'/'pca'/'mds'
do_h_cluster      = 'ye';               % 'yes'/'no' Hierarchical Clustering ?
do_kmeans_cluster = 'ye';              % 'yes'/'no' K-means clustering?
% --------------> Data
data_mode         =  'feature_vector';  % 'feature_vector'/'distance_matrix'
data_file         = '.\human_orig.txt';         % feat., complete
                    % -> set 0/999 to use all numerical data (starting at cell (2,2)) ->
  col_varfirst    =  0;                 % 0    to start with the first numerical column
  col_varlast     =  9999;              % 9999 or any large number to use til last column
norm_mode         =  2;                 % normalize 0/ none       1/divide by SD
                                        %           2/z-scores;   3/div. by samp. mean
                                        %           4/grand-means 5/grand-medians
% --------------> Hierarchical cluster analysis
save_cluster    = 'no';                 % Save cluster number of each case?        ['no']
dist_metric     = 'euclid';             %>Also used in MDS/K-Means-clust. Features only
                                        %   'euclid'/'seuclid'/'cosine'/'correlation'/
                                        %   'cityblock'                        ['euclid']
                                        %>Distance-matrix only: select 'euclid'/'seuclid'
                                        % Must conform the type of inp. distance matrix
link_method     = 'ward';               % 'single'/'complete'/'average'/'centroid'/'ward'
% --------------> Multidimensional scaling (MDS)
                  % The distance matric set in /Hierarchical cluster analysis/ is used
mds_criterion     = 'stress';           % 'stress'/'metricstress'              ['stress']
% --------------> Graphic layout
maxn_numbers      =  60;                % maximum number of data points to be labeled
FONT_SIZE         =  8;
%>}>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% not! to be changed by user!
%data_file        = '.\test_data\digrju_dist_mat_euc1_cluster.txt';   % Distance matrix
%                    Path-name to grand-means reference data is fixed to ->
grand_means       =  [MLMNETPATH '\data\grand_mean_cortex.csv'];
%  cutoff_level      Resolution of clusters, relative cut-off level  (0..1)
%                    0.7/ default; 0.0/ default; 1/ only one cluster
%                    Will be replaced by an estimated value (spss) ->
cutoff_level      =  0.7;          % Relative cut-off level  (0.0 ... 1.0)          [0.7]
do_transpose      = 'no';          % 'yes'/'no'  not yet implemented
activate_stop     = 'no';          % 'ye'/'no' activates <message.window> to show
                                   %           figures window
i_figure   = 0;
break_line(1:80)  = '-';
%--------------------------------------------------------------------------------- m+read
  [err,raw_data,var_names,case_names] = ...
    y_tblread_block(data_file,col_varfirst,col_varlast);
% data transform (optional, to be extended)
  if 0; raw_data = sqrt(raw_data); end

  [m,n]=size(raw_data);
  N_ROWS    = m; N_COLUMNS = n;
  number_variables = n;                     % use all columns!
  s_text = ['\n.... ASCII-File <%s> loaded: \n' ...
             '     (%3.0f) rows, (%3.0f) columns of numerical data\n'];
  fprintf(1,s_text,data_file,N_ROWS,N_COLUMNS);
  case_names_c =cellstr(case_names);   % creat CAOS for case(=line) names
  [case_names_number,case_names_length] = size(case_names);

  % Save Input data as *.mat file (and reload)
    [pathstr,file_name,file_ext] = fileparts(data_file);
    file_out = [pwd '\OUT\' file_name '.mat'];
    save(file_out,'raw_data','var_names','case_names'); % load(file_out);
  % Check if input data is a square distance matrix (not trimmed!)
    if strncmpi(data_mode,'di',2)
     if ( N_ROWS == N_COLUMNS) && ...
            (raw_data(N_ROWS,1) == raw_data(1,N_COLUMNS))
             fprintf(1,'.... Inputfile accepted as %3.0f x%3.0f distance matrix\n', ...
             N_ROWS,N_COLUMNS);
     else
             error('???? Inputfile is not a square distance matrix'); beep; beep
     end
    end

% { ----------------------------------------------------------------------------- m+norm
%  Normalization of raw data, not to be used with distance matrices!
%  Normalized data  will be stored to <.\OUT using tbl-format>
   if strncmpi(data_mode,'fe',2)    %----{ Normalization, features only

       % Normalize by grand-mean/median (optional)                         % m+grand-mean
       if ((norm_mode == 4) || (norm_mode == 5)) %----{ Normalize by 4/grand-mean 5/median
       % Read data
       % Results are
       % * lgmvar_names  (lower case, length reduced to 4 characters)
       % *  gm_data      Columns of grand-mean data file which correspounds to input file
         [err,gmraw_data,gmvar_names,gmcase_names] = ...
             y_tblread_block(grand_means,0,9999);
       % Check data. Variable names set to lowercase. Only the first 4 chars are used!
       % Feature data: variable names [ligand names]  <var_names>,  Char array of strings
       % Grand means:  variable names [ligand names]  <gmvar_names>,Char array of strings
         type(grand_means); % show complete grand-means data file
         lvar_names   = lower(var_names(:,1:4));     % var-names of data to be analyzed
         lgmvar_names = lower(gmvar_names(:,1:4));   % var-names of grand-means data file
         compare_names = strcmp(lvar_names, lgmvar_names);
         if (compare_names)
           % Number, content and sequence of variable names are identical
           % No action required
         else
           % Number, content and sequence of variable names differ
           % Action: Pick up the corresponding columns indicated in <var_names>
           %         from 1) gmvar_names and 2) gmraw_data
           [m_data,n_data] = size(lvar_names);    % data to be processed, n_data
           [m_gm,  n_gm]   = size(lgmvar_names);  % reference data
           if (m_data <= m_gm)
               % select and reorder data from <gm: grand mean/med> to fit input data
               entries_found = zeros(m_data,1);   % col vector
               for i_run = 1:m_data
                   try
                    entries_found(i_run)=strmatch(lvar_names(i_run,:),lgmvar_names);
                   catch
                    error(['    Variable <' lvar_names(i_run,:) ...
                           '> not found in the gm-file']);
                   end
               end
               if length(entries_found) < m_data
                  error(' Some variable names not found in the gm-file')
               end
               gmraw_data  = gmraw_data(:,entries_found);
               gmvar_names = gmvar_names(entries_found,:);
           else
               % more variables in input file than in the reference file (grand_mean)
               % stop procedure
                error(' More variable names than in the gm-file')
           end
         end
       end                                      %----}     Normalize by grand-mean/median

       switch norm_mode% assign descriptive names to normalization modes, no calculations
        case 0,    snorm_mode  = 'none';                     name_app = '_none';
        case 1,    snorm_mode  = 'divide-by-SD';             name_app = '_dbSD';
        case 2,    snorm_mode  = 'z-scores';                 name_app = '_zsco';
        case 3,    snorm_mode  = 'divide-by-sample mean';    name_app = '_dbsmean';
        case 4,    snorm_mode  = 'divide-by-grand mean';     name_app = '_dbgmean';
                   grand_mean  = gmraw_data(1,:);
        case 5,    snorm_mode  = 'divide-by-grand median';   name_app = '_dbgmedian';
                   grand_median= gmraw_data(2,:);
        otherwise, snorm_mode  = 'none';                     name_app = '_none';
       end

       switch norm_mode
         case 0 %- No normalization, use raw data
                   norm_data = raw_data;     % just copy data
                   stand_text = 'No standardization';
         case 1 %- Divided by SD)          Normalize
                   [m,n] = size(raw_data);
                   stdr = std(raw_data);
                   norm_data = raw_data./stdr(ones(m,1),:);
                   stand_text = 'Normalization: Divide by SD';
                   normal_colstat = mean(norm_data);
         case 2 %- z-scores,               Normalize
                   norm_data = zscore(raw_data);
                   stand_text = 'Normalization: z-Scores';
                   normal_colstat = mean(norm_data);
         case 3 %- Divided by sample mean  Normalize
                   [m,n] = size(raw_data); mean_values = mean(raw_data);
                   norm_data = raw_data./mean_values(ones(m,1),:);
                   stand_text = 'Normalization: Divide by sample mean';
                   normal_colstat = mean(norm_data);
         case 4 %- Divided by grand mean   Normalize
                   [m,n] = size(raw_data);
                   norm_data = raw_data./grand_mean(ones(m,1),:);
                   mgrand_means = strrep(grand_means,'\','/');
                   stand_text = 'Normalization: Divide by grand mean';
                   % Calculate mean of normalized variables.
                   normal_colstat = mean(norm_data);
         case 5 %- Divided by grand median Normalize
                   [m,n] = size(raw_data);
                   norm_data = raw_data./grand_median(ones(m,1),:);
                   mgrand_means = strrep(grand_means,'\','/');
                   stand_text = 'Normalization: Divide by grand median';
                   % Calculate median of normalized variables.
                   normal_colstat = mean(norm_data);
         otherwise %- Number not valid, use raw data, no standardization
                   norm_data = raw_data;
                   stand_text = 'No standardization, feature vectors';
       end      % end switch
   else
       norm_data = raw_data;
       stand_text = 'No standardization, data is a distance matrix';
   end                               %----} Normalization, features only

 % Output of  data (raw or normalized), tbl-format using case- and variable-names
   s_text = ['.... ' stand_text '\n']; fprintf(1,s_text);
   [pathstr,file_name,file_ext] = fileparts(data_file);
   file_out = [pwd '\OUT\' file_name name_app '.csv'];
   tblwrite(norm_data,var_names,case_names,file_out,'space')
   % Information on saving standardized data -> screen:
   fprintf(1,'     (Normalized) data saved:\n     %s\n',file_out);
   format short; type(file_out)

  % Output of variable means or medians to indicate the variable specific weight
  % tbl-format using variable-names
    if (norm_mode >0)
       % set values of <normal_colstat> near 1 to zero
       zero_index = find(abs(normal_colstat) <0.1e-8);
       normal_colstat(zero_index) = 0 ;
       [pathstr,file_name,file_ext] = fileparts(data_file);
       file_out = [pwd '\OUT\' file_name name_app '_colstat.csv'];
       tblwrite(normal_colstat,var_names,'colstat',file_out,'space')
       % Information on saving standardized data -> screen:
       fprintf(1,'     Column statistics (means or medians) saved:\n     %s\n',file_out);
       type(file_out)
    end

% } --------------------------------------------------------------------------------------

% { -------------------------------------------------------------------------- m+transpose
if strncmpi(data_mode,'fe',2)
   if (strncmpi(do_transpose,'ye',2)) % transpose, result: rows=features, cols=objects
   norm_data = (norm_data)';
   fprintf(1,'\n.... data transposed!\n');
   xlabel_text = 'Value';
   ylabel_text = 'Variables (columns)';
   else % use original data; rows = objects, cols=features
      % fprintf(1,'.... Data not transposed\n');
      xlabel_text = 'Variables (columns)';  ylabel_text = 'Value';
   end
end
% } --------------------------------------------------------------------------------------

% { ---------------------------------------------------------------------------- m+boxplot
% Boxplot(X) produces a box and whisker plot for each column of X.
% For feature vectors only!
%   The box has lines at the lower quartile, median, and upper quartile values.
%   The whiskers are lines extending from each end of the box to show the extent of the
%   rest of the data. Outliers are data with values beyond the ends of the whiskers.
%   If there is no data outside the whisker, a dot is placed at the bottom whisker.

if (strncmpi(do_boxplot,'ye',2)) && (strncmpi(data_mode,'fe',2))
  i_figure = i_figure+1;
  figure(i_figure); clf; z_setgraph
  boxplot(norm_data,'labels',var_names,'orientation','vertical')
  set(gcf,'DefaulttextInterpreter','none')
  [pathstr,name,ext] = fileparts(data_file);
  xlabel([xlabel_text '  /  ' stand_text]);
  ylabel(ylabel_text);
  my_text = sprintf('%s%s / %s / Transpose:',name,ext,datestr(now),do_transpose);
  my_text_new = strrep(my_text,'_','-'); title(my_text_new);
  if 0
    if (strncmpi(do_transpose,'ye',2)); set(gca,'YTicklabel',case_names)
    else                                set(gca,'YTicklabel',var_names)
    end
  end
  y_save_fig(PATH_OUT,file_name,'_boxplot',GRAPH_FORMAT,[0;0]);
  if strncmpi(activate_stop,'ye',2)
      beep; k=y_showgraph;
      if (k == 1); print -dwinc; end
  end
 end
% } --------------------------------------------------------------------------------------

% { -------------------------------------------------------------------------------- m+pca
%>>Calculate  principal components and generate plot: princomp()
%  For feature vectors only!
%  pcs:       principal components: new axes in old system
%  pca_data:  datapoints in the new coordinate system
%  variances: variance explained by the corresponding column
%             (principal component)
%  t2         Distance of datapoints from center
%  Problem:   Number of variables cannot exceed the number of rows (cases)
%             Use pcacov() to calculate scores from covariance matrix
%  INPUT:     norm_data, no distance matrix accepted
%  OUTPUT:    pca_scores
%----------------------------------------------------------------------------------------
if (strncmpi(do_scatter_plot,'pc',2)) && strncmpi(data_mode,'fe',2)
   covx = cov(norm_data);                   % Covariance matrix
   [pc,variances,explained] = pcacov(covx); % PCA from Covariance matrix

   % print the first four principal components
     [m,n_pca]=size(pc); if (n_pca > 4);n_pca = 4; end
     s_text = ['\n.... PCA: First %3.0f principal components\n' break_line '\n'];
     fprintf(1,s_text,n_pca);
     first_pc = pc(:,1:n_pca); [m,n] = size(first_pc);
     for i_run = 1:m;
        fprintf(1,'\t%8.4f %8.4f %8.4f %8.4f -- %s \n', ...
                   (first_pc(i_run,:))',var_names(i_run,:))
     end

   cum_explained    = cumsum(explained); 
   cum_explained_2d = cum_explained(2);     % varaince explained in a 2d-graph
   s_text = ['\n.... PCA: percentage of the total variance explained by each '...
                         'eigenvector: \n' break_line '\n'];
   running_number = (linspace(1,number_variables,number_variables))';
   fprintf(1,s_text);
   fprintf(1,'      %2.0f %8.3f  %8.3f\n', [running_number explained cum_explained]')

   col_means = mean(norm_data);             % means of variables
   [m,n] = size(norm_data);
   col_means_mat = col_means(ones(m,1),:);  % Expand col-vector to matrix
                                            % Shift data to origin=translation
   norm_data_transl = (norm_data-col_means_mat);
   pca_scores = (norm_data_transl)*pc;      % PCA

   s_text = ['\n.... PCA-scores calculated. First values are: \n' break_line '\n'];
   fprintf(1,s_text); [m,n] = size(pca_scores);
   up_to = 4; if (up_to) > n; up_to = n; end
   pca_scorest = pca_scores(1:4,1:up_to)';
   fprintf(1,'      %8.4f %8.4f %8.4f %8.4f \n',pca_scorest ); fprintf(1,'\n');
   first_pc  = pca_scores(:,1); second_pc = pca_scores(:,2); i_figure = i_figure+1;
   % check range of axes
     grand_min = min(min(pca_scores(:,1:2)));
     grand_max = max(max(pca_scores(:,1:2)));
     x_shift = (grand_max-grand_min)/100;
   % Make PCA-plot
     figure(i_figure); clf; z_setgraph
     plot (first_pc,second_pc,'b+',first_pc,second_pc,'ro')
     set(gcf,'DefaulttextInterpreter','none')
     xlabel('1st Principal Component'); ylabel('2nd Principal Component')
     title_line_1 = ['=PCA= ' data_file ' / ' datestr(now) '  ' stand_text];
     title_line_2 = ['Variance explained in 2D: ' num2str(cum_explained_2d,'%8.3f') '%'];
     title_array =  {title_line_1;title_line_2};
     title (title_array);
     axis square; set(gca,'FontSize',FONT_SIZE)
     %grid on
   % Label datapoints by row-number, number limited by <maxn_numbers>
     if (N_ROWS <= maxn_numbers)
        for i = 1:m
              s_text = case_names(i,:);
              text(first_pc(i)+x_shift,second_pc(i)-x_shift,s_text,'FontSize',FONT_SIZE)
        end
     end
     beep; k=y_showgraph;
   % save graph PCA as *.emf file
     y_save_fig(PATH_OUT,file_name,'_pca',GRAPH_FORMAT,[0;0]);
     if strncmpi(activate_stop,'ye',2)
      beep; k=y_showgraph;
      if (k == 1); print -dwinc; end
     end

end
% } -------------------------------------------------------------------------------- m+pca

% { -------------------------------------------------------------------------------- m+mds
%>> Multidimensional scaling. Inputdata:  norm_data
%   Each line (row) in <norm_data> is one case (object)
%   A valid value for <dist_metric> must be defined (see m+'user').
%   Data with distance_matrix accepted
%   ----------------------------------------------------------------
%>> -cmdscale() performs classical multidimensional scaling (MDS).
%   -mdscale()  performs nonclassical MDS. As with cmdscale, use mdscale either to
%   visualize dissimilarity data for which no "locations" exist, or to visualize
%   high-dimensional data by reducing its dimensionality. Both take a matrix of
%   dissimilarities as an input and produce a configuration of points. mdscale offers a
%   choice of different criteria to construct the configuration, and allows missing data
%   and weights. -> mdscale(): Nonmetric and metric multidimensional scaling
%   -  By default, mdscale() uses Kruskal's normalized stress1 criterion
%   -  'Criterion'?The goodness-of-fit criterion to minimize:
%      'stress' ?Stress normalized by the sum of squares of the interpoint distances,
%                 stress1. Default. Non-metric scaling !!
%      'metricstress' ?Stress normalized with the sum of squares of the
%                 dissimilarities.  Metric scaling !!
if strncmpi(do_scatter_plot,'md',2)
   s_text1 = ...
        '\n.... MDS (Multidimensional scaling).  Metric: %s.  Criterion: %s\n';
   fprintf(1,s_text1,dist_metric,mds_criterion);
   if strncmpi(data_mode,'fe',2)
      % Input are feature vectors, calculate distances
        dissimilarities   =  pdist(norm_data,dist_metric); % dissimilarities: a vector!
   else
      % Input is a distance matrix. From square matrix to vector
        dissimilarities = squareform(norm_data);
   end
   %mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
   opts = statset('MaxIter',2000');   %  Max number of iterations. Default = 200
   [points_2d,stress,disparities] = ...
        mdscale(dissimilarities,2,'Criterion',mds_criterion,'Options',opts);
   %mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
   fprintf(1,'     Stress: %6.4f\n',stress);

% Shepard plot of the <distances>.
  distances = pdist(points_2d,dist_metric);  % distances between points in 2-d
  [dum,ord] = sortrows([disparities(:) dissimilarities(:)]);
  i_figure = i_figure+1;
  figure(i_figure); clf; z_setgraph
  hold on  % ppppppppppppppppppppppppppppppppppppppppppppppppp
      plot(dissimilarities,distances,'bo', ...
               dissimilarities(ord),disparities(ord),'r-','MarkerSize',4);
      set(gcf,'DefaulttextInterpreter','none')
      xlabel('Dissimilarities [observed]');
      ylabel('Distances [estim.] / Disparities [red]')
      legend({'Distances' 'Disparities'}, 'Location','NorthWest');
      title  (['=MDS= Shepard plot' data_file ' / ' datestr(now) ' / ' dist_metric ...
                      ' / ' mds_criterion ': ' num2str(stress,'%8.6f')])
        [x_limits,y_limits]=y_get_axes_limits(gcf); % plot diagonal
        plot (x_limits',y_limits','k-');
  hold off %pppppppppppppppppppppppppppppppppppppppppppppppppp

% Make MDS-Plot
  i_figure = i_figure+1; figure(i_figure); clf; z_setgraph
   grand_min = min(min(points_2d(:,1:2)));
   grand_max = max(max(points_2d(:,1:2)));
   x_shift = (grand_max-grand_min)/100;
   plot(points_2d(:,1),points_2d(:,2),'b+', ...
        points_2d(:,1),points_2d(:,2),'ro');
   set(gcf,'DefaulttextInterpreter','none')
   set(gca,'FontSize',FONT_SIZE); set(gcf,'DefaulttextInterpreter','none')
   case_names_c =cellstr(case_names);
   text(points_2d(:,1)+x_shift,points_2d(:,2)-x_shift,case_names_c, ...
        'FontSize',FONT_SIZE);
   title  (['=MDS='  data_file '  /  ' datestr(now) '  /  ' dist_metric ...
           ' / ' mds_criterion ': ' num2str(stress,'%8.6f')])
   xlabel('Dimension 1'); ylabel('Dimension 2'); %'
   axis square; axis on
   % equal scale for both axes
     [x_limits,y_limits]=y_get_axes_limits(gcf);
     lower_xy=min([x_limits(1) y_limits(1)]);
     upper_xy=max([x_limits(2) y_limits(2)]);
     xlim([lower_xy upper_xy]); ylim([lower_xy upper_xy]);
   % save graph MDS as *.emf file
     y_save_fig(PATH_OUT,file_name,'_mds',GRAPH_FORMAT,[0;0]);
     if strncmpi(activate_stop,'ye',2)
      beep; k=y_showgraph;
      if (k == 1); print -dwinc; end
     end

end % END: if strncmpi(scatter_plot,'md',2)
% } -------------------------------------------------------------------------------- m+mds

% { ---------------------------------------------------------------------------- m+cluster
%>> Hierarchical cluster analysis. Data processed: norm_data                              #+mm
%   Distance_matrix input also accepted
if strncmpi(do_h_cluster,'ye',2)

    if strncmpi(data_mode,'fe',2)  
	    % input data are 'features'
         distances    = pdist(norm_data,dist_metric);  % vector                  #+m-pdist
         %distances   = pdist(pca_scores(:,1:end),dist_metric); % Results identical!
         dist_matrix  = squareform(distances);         % matrix
         % output distance matrix. Dist-metric given by str(dist_metric,1,3)
         [pathstr,file_name,file_ext] = fileparts(data_file); % data input
         file_out = [pwd '\OUT\' file_name '_' dist_metric(1:3) '-dist.csv'];
         tblwrite(dist_matrix,case_names,case_names,file_out,'tab');
         fprintf(1,'.... Distance-matrix saved:\n     %s\n',file_out);
         format short e; type(file_out)  % writes dist-matrix to screen
    else                           
	     % input data is a 'distance matrix'
         distances    = squareform(norm_data);          % vector
         dist_matrix  =            norm_data;           % matrix
    end
	
    n_dendro = 0;
    dist_matrix,        % Display

while 1  % { START >  while (modify dendrogram)
   n_dendro = n_dendro+1;                       % counts dendrograms
   mk_clusters = linkage(distances,link_method); %links between objects:       #+m-linkage
  % Each dendrogram() creates its own figure window. All elements are plotted
  % Distances between fused objects in <mk_clusters(:,3)

    if n_dendro == 1;      % { start initial cluster number, first dendrogram only ----
      % Analyse sequence of distances. Estimates an initial cluster number (see spss)
      merge_diff   = [diff(mk_clusters(:,3));0];
      mk_clusters  = [mk_clusters merge_diff];
      max_diff     = max(merge_diff); ind_max_diff = find(merge_diff == max_diff);
      dist_level   = mk_clusters(ind_max_diff:ind_max_diff+1,3);
      ini_cutoff   = (dist_level(1)+(dist_level(2)-dist_level(1))/2)/ ...
                      max(mk_clusters(:,3));
      cutoff_level = ini_cutoff;
    end                    % } end initial number of clusters ----

  % Plot dendrogram
    i_figure = i_figure+1; figure(i_figure); clf; z_setgraph
    % Setting the cutoff-level' ( -> number of clusters)
      max_mk_clusters_3 = max(mk_clusters(:,3));        % distance highest level
      color_threshold = cutoff_level*max_mk_clusters_3; % an absolute distance
      %mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
      [H,T,perm] = dendrogram(mk_clusters,0,'colorthreshold',color_threshold, ...
                              'orientation','right');
  %-Modify graph, label y-axis
    set(gcf,'DefaulttextInterpreter','none')
    corr_coefficient = cophenet(mk_clusters,distances);
    s_text1=sprintf('%s  /  %s',datestr(now),data_file);
    s_text2=sprintf('%s  /  %s  / %s  / Cophenetic: %6.4f', ...
                    stand_text,dist_metric,link_method,corr_coefficient);
    s_array = {s_text1,s_text2};
    title (s_array);  xlabel('Distance'); ylabel('Cases');
    n_ticks = length(case_names_c);
    for i_tick = 1:n_ticks
       tick_label{i_tick}=case_names_c{perm(i_tick)};
    end
    set(gca,'YTickLabel',tick_label)
    set(gca,'Box','On');  set(gca,'FontSize',FONT_SIZE)
    axis square; axis on ;

    if n_dendro == 1;        % start first dendrogram only ----
        % output only for the first dendrogram
        %-graph title
          s_text = [  '.... Hierarchical cluster analysis using:' ...
                      '\n     Metric: ' dist_metric ', linkage: ' link_method ...
                      '\n     Cophenetic correlation coefficients is :%8.5f'  ...
                      '\n'];
          fprintf(1,s_text, corr_coefficient);
        %-Linkage
          if 0 %debug
             s_text = '.... Linkage information /object-1/object-2/distance/:\n';
             fprintf(1,s_text);
             fprintf(1,'    %4.0f %4.0f %9.4f\n',(mk_clusters)');
          end
        %-Inconsistency coefficients. maximum values only
          inconsist_mat = inconsistent(mk_clusters);
          % n_inconsistencies = 5;
          % inconsist_mat_s = sortrows(inconsist_mat,4);
          % s_text = '     Max inconsistencies: /Mean height/SD/Number/Coefficient/:\n';
          % fprintf(1,s_text);
          % [m,n] = size(inconsist_mat_s);
          % if m < n_inconsistencies;   % print all
          %    fprintf(1,'    %9.4f %9.4f %4.0f %9.4f\n',(inconsist_mat_s(:,:))');
          % else
          %    fprintf(1,'    %9.4f %9.4f %4.0f %9.4f\n', ...
          %              (inconsist_mat_s(end-n_inconsistencies:end,:))');
          % end
    end                     % end first dendrogram only ----

  % indicate active cutoff_level and number of clusters
    dummy_1 = round(cutoff_level*max_mk_clusters_3*100)/100;            % absolute
    %cutoff_level     = str2num(c_array{1})/ max_mk_clusters_3;         % relative 0..1
    cutoff_level_abs  = round(cutoff_level*max_mk_clusters_3*100)/100;  % absolute
    fprintf(1,'   --Cutoff level (rel/abs) set to %6.4f / %6.4f\n', ...
              cutoff_level,cutoff_level_abs);
    % get detected number of clusters;
    n_hcluster = length(find(mk_clusters(:,3) > cutoff_level_abs))+1;
    fprintf(1,'     Number of clusters detected: %5.0f\n',n_hcluster);
  %-get cluster number of each case
    cluster_index = cluster(mk_clusters(:,1:3),'maxclust',n_hcluster);
  % Discriminant analysis using classify(), to validate results of H-clustering
    try
        [class,err,posterior] = classify(norm_data,norm_data,cluster_index,'linear');
        %class, %posterior
        s_text = [ '     * Misclassification error rate using\n' ...
                   '       Classify(), linear discriminant function (%%): %6.4f \n'];
        fprintf(1,s_text,err*100);
    catch
    end
  %-Enter new cutoff_level -> <cutoff_level>, absolute distance
    prompt = 'New cutoff level (distance), [Cancel] to continue';
    dlg_title  = 'Cutoff level'; beep;
    de_fault = {num2str( dummy_1)};
    c_array = inputdlg(prompt,dlg_title,1,de_fault);
    if (isempty(c_array) || str2num(c_array{1}) == 0); break; end       % Cancel skips
    cutoff_level = str2num(c_array{1})/ max_mk_clusters_3;              % relativ
end      % } ENDE  <  while (modify dendrogram)

  % interrupt to show graph, not if last plot
    if strncmpi(activate_stop,'ye',2)
      beep; k=y_showgraph;
      if (k == 1)
	     print -dwinc
	  end
    end
	
  % save graph Cluster as *.emf file
    y_save_fig(PATH_OUT,file_name,'_hclus',GRAPH_FORMAT,[0;0]);

  % Save cluster number of each case
    if strncmpi(save_cluster,'ye',2)
      % save cluster-indices
        [pathstr, name, ext] = fileparts(data_file); % get name of input file
        % cases-names: <case_names>; cluster-indices: <cluster_index>
        file_out =['.\' name '_clusterind.txt'];
        [m,n] = size(cluster_index); fid = fopen(file_out,'w+');
        for i_run =1:m
		    fprintf(fid,'%-16s %5.0f\n',case_names(i_run,:),cluster_index(i_run));
		end
        fprintf(1,'.... H-Cluster: Cluster indices saved to file %s\n',file_out);
    end
	
end        % end CLUSTER
% } ---------------------------------------------------------------------------- m+cluster

% { ---------------------------------------------------------------------------- m+k-means
%>>K-means cluster analysis. Inputdata: norm_data
%  Only for 'feature_vector' data
if strncmpi(do_kmeans_cluster,'ye',2) && strncmpi(data_mode,'fe',2)   %--{
   % Check for valid distance measures
           valid_dist_meas = {'none';'sqEuclidean';'cityblock';'cosine';'correlation'};
           if strcmpi(dist_metric,'euclid'); dist_metric = 'sqEuclidean'; end
           entry_found = strmatch(dist_metric,valid_dist_meas);

   if (entry_found)
     n_replicates = 8;
     % hierarchical- and 'k-means' clustering use different names for Euclidean (squared)
       n_cluster = 2; c_array = {'2'};
       i_figure = i_figure+1;
       figure(i_figure); clf; z_setgraph
	   
    while 1    % --{
       prompt =  ...
       'Enter a valid number number of clusters (> 1). [Cancel] to continue';
       dlg_title  = 'K-means clustering'; beep;
       de_fault = {int2str(n_cluster)};
       c_array = inputdlg(prompt,dlg_title,1,de_fault);
       if (isempty(c_array) | str2num(c_array{1}) < 2); break; end   % Cancel skips
       n_cluster = str2num(c_array{1});
       s_text = '.... K-means clustering,%3.0f clusters,%3.0f replicates:\n';
       fprintf(1,s_text,n_cluster,n_replicates);
       % kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
       [idx_kmeans,clust_cent,clust_sumdis,point_distances]= ...
           kmeans( ...
           norm_data,n_cluster,'dist',dist_metric,'display','final', ...
           'replicates',n_replicates);
       % Silhouette-plot
         [silh_ncluster,h] = silhouette(norm_data,idx_kmeans,dist_metric);
         mean_silh =mean(silh_ncluster);
         set(gcf,'DefaulttextInterpreter','none')
         s_text=sprintf('%s  /  %s / %s  /  %s  / Separation:%6.4f', ...
                        data_file,datestr(now),stand_text,dist_metric, ...
                         mean_silh);
         title (s_text)
       % Silhouette-Mean
         s_text = '\n     Mean silhouette value using  %3.0f clusters: %6.5f\n';
         fprintf(1,s_text,n_cluster,mean_silh);
       % Show cluster members
	   
        for idx_cluster = 1:n_cluster     % idx_cluster: Loop over clusters------------
              fprintf(1,'%s\n',s_line);
              % <this_cluster> : indices of all members of this cluster
              this_cluster = find(idx_kmeans == idx_cluster);
              n_members=length(this_cluster);
              fprintf(1,['     Cluster %2.0f contains %4.0f cases' ....
                         ' (Distance to center/Case name)\n'],idx_cluster,n_members);
              % <this_cluster_dist>: dist. of all members of this cluster to centroid
              this_cluster_dist = point_distances(this_cluster,idx_cluster);
              % i_this_cluster: Loop over members in cluster---
			  
            for i_this_cluster = 1:n_members
                fprintf(1,'    %8.4f   %s   \n', ...
                          this_cluster_dist(i_this_cluster)  , ...
                          case_names_c{this_cluster(i_this_cluster)});
            end
			  
        end
		  
        fprintf(1,'\n')
        % save silhouette-cluster graph (active figure)
        y_save_fig(PATH_OUT,file_name,'_kclus',GRAPH_FORMAT,[0;0]);
			  
    end % end WHILE   % -- }
		  
        if strncmpi(activate_stop,'ye',2)
                beep; k=y_showgraph;
                if (k == 1); print -dwinc; end
        end
         
        fprintf(1,'???? <%s> is not a valid dis- measure for a cluster analysis\n', ...
                    dist_metric)
     end     % end k-means
end         % end -> if strncmpi(do_kmeans_cluster,'ye')  %--}
% } ---------------------------------------------------------------------------- m+k-means

diary off   % Close diary
% } --------------------------------------------------------------------------- m+main end

function  []=z_setgraph()
%--------------------------------------------------------------------------- m+z_setgraph
% INPUT ARGUMENTS:  none
% OUTPUT ARGUMENT: none
%
% Definition of Plotsize (see Manual: Using Graphics, p9-20) of the current figure:
%----------------------------------------------------------------------------------------
%  v(1) Lower left corner: distance to right margin
%  v(1) Lower left corner: distance to bottom
%  v(1) height of plot:    (in x)
%  v(1) width  of plot:    (in y)
%----------------------------------------------------------------------------------------
set(gcf,'Paperorientation', 'portrait')      % Defines paper orientation
set(gcf,'Paperunits', 'centimeters')         % Defines scale units
set(gcf,'PaperType',  'A4')                  % Defines Paper size
set(gcf,'Units','pixels')
%- Define Plot-Window:
   screen_size = get(0,'Screensize');
   pos_vec = [screen_size(3)/5    screen_size(4)/5 ...
                          screen_size(3)*0.75 screen_size(4)*0.70];
   set(gcf,'Position',pos_vec)  %set size in current figure
   h=axes('Position',[0.15 0.1 0.80 0.80]);  % defines borders
   set(gcf,'Interruptible','on')
%- Define Harcopy features for prints
   print_pos = [2.0 7.0 18.0 18.0];
   set(gcf,'Paperposition',print_pos)
