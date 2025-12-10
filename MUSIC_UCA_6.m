classdef MUSIC_UCA_6 < handle
    % MUSIC_UCA_6 - 6-Element Uniform Circular Array MUSIC Algorithm
    % ===============================================================
    %
    % Optimized MUSIC AOA estimation for 6 bowtie dipole antennas
    % with parabolic reflectors in UCA configuration.
    %
    % System Specifications:
    % - Array: 6 elements, radius 0.176m (0.6λ at 1050 MHz)
    % - Frequency Range: 700-1400 MHz
    % - Directional antennas: 60° HPBW, 15 dB F/B ratio
    %
    % Author: Ported from Python implementation
    % Date: 2025-12-08

    properties
        % Array configuration
        num_elements = 6
        radius_m = 0.176  % meters (0.6λ at 1050 MHz)
        center_freq_hz = 1050e6  % Hz

        % Antenna pattern (wider beamwidth for better coverage)
        antenna_beamwidth_deg = 120.0  % HPBW (wider for 6-element UCA)
        antenna_fb_ratio_db = 15.0  % Front-to-back ratio (increased to suppress ghost peaks)

        % Algorithm parameters (increased to 4 for ghost/multipath handling)
        num_sources = 4

        % Computed properties
        antenna_positions  % (6, 3) array [x, y, z] in meters
        antenna_orientations  % (6,) azimuth pointing directions
        wavelength_m
    end

    methods
        function obj = MUSIC_UCA_6(varargin)
            % Constructor
            %
            % Usage:
            %   music = MUSIC_UCA_6()
            %   music = MUSIC_UCA_6('num_sources', 3, 'frequency', 1050e6)

            % Parse optional inputs
            p = inputParser;
            addParameter(p, 'num_sources', 3);
            addParameter(p, 'frequency', 1050e6);
            addParameter(p, 'radius', 0.176);
            parse(p, varargin{:});

            obj.num_sources = p.Results.num_sources;
            obj.center_freq_hz = p.Results.frequency;
            obj.radius_m = p.Results.radius;

            % Compute wavelength
            c = 3e8; % speed of light
            obj.wavelength_m = c / obj.center_freq_hz;

            % Initialize antenna geometry
            obj.compute_antenna_geometry();

            fprintf('MUSIC_UCA_6 initialized:\n');
            fprintf('  Elements: %d\n', obj.num_elements);
            fprintf('  Radius: %.4f m (%.2fλ)\n', obj.radius_m, ...
                obj.radius_m / obj.wavelength_m);
            fprintf('  Frequency: %.1f MHz\n', obj.center_freq_hz / 1e6);
            fprintf('  Expected sources: %d\n', obj.num_sources);
        end

        function compute_antenna_geometry(obj)
            % Compute antenna positions and orientations

            % Antenna positions (planar UCA in X-Y plane)
            angles = linspace(0, 2*pi, obj.num_elements + 1);
            angles = angles(1:obj.num_elements);

            obj.antenna_positions = zeros(obj.num_elements, 3);
            obj.antenna_positions(:, 1) = obj.radius_m * cos(angles);  % x
            obj.antenna_positions(:, 2) = obj.radius_m * sin(angles);  % y
            obj.antenna_positions(:, 3) = 0;  % z (planar)

            % Antenna orientations (pointing toward center)
            obj.antenna_orientations = mod(rad2deg(angles) + 180, 360);
        end

        function gain = directional_antenna_gain(obj, azimuth_deg, antenna_idx)
            % Calculate directional antenna gain (Gaussian beam pattern)
            %
            % Inputs:
            %   azimuth_deg: Signal DOA azimuth (degrees)
            %   antenna_idx: Antenna index (1-6)
            %
            % Returns:
            %   gain: Linear scale gain (0-1)

            % Angle difference between signal and antenna pointing
            antenna_pointing = obj.antenna_orientations(antenna_idx);
            angle_diff = abs(azimuth_deg - antenna_pointing);

            % Wrap to [-180, 180]
            if angle_diff > 180
                angle_diff = 360 - angle_diff;
            end

            % Gaussian beam pattern
            sigma = obj.antenna_beamwidth_deg / (2 * sqrt(2 * log(2)));

            % Apply front-to-back ratio
            if abs(angle_diff) > 90  % Back lobe
                gain = 10^(-obj.antenna_fb_ratio_db / 20) * ...
                    exp(-(angle_diff - 180)^2 / (2 * sigma^2));
            else  % Main lobe
                gain = exp(-angle_diff^2 / (2 * sigma^2));
            end
        end

        function a = steering_vector(obj, azimuth_deg, elevation_deg, ...
                include_antenna_pattern)
            % Calculate UCA steering vector
            %
            % Inputs:
            %   azimuth_deg: Azimuth angle (degrees, -180 to 180)
            %   elevation_deg: Elevation angle (degrees, default 0)
            %   include_antenna_pattern: Include directional pattern (default true)
            %
            % Returns:
            %   a: (6×1) complex steering vector

            if nargin < 3
                elevation_deg = 0;
            end
            if nargin < 4
                include_antenna_pattern = true;
            end

            % Convert to radians
            az_rad = deg2rad(azimuth_deg);
            el_rad = deg2rad(elevation_deg);

            % Wave number
            k = 2 * pi / obj.wavelength_m;

            % Direction of arrival vector (unit vector)
            doa = [
                cos(el_rad) * cos(az_rad);
                cos(el_rad) * sin(az_rad);
                sin(el_rad)
            ];

            % Calculate phase shifts
            phase_shifts = k * (obj.antenna_positions * doa);

            % Base steering vector
            a = exp(1j * phase_shifts);

            % Apply directional antenna gain
            if include_antenna_pattern
                antenna_gains = zeros(obj.num_elements, 1);
                for i = 1:obj.num_elements
                    antenna_gains(i) = obj.directional_antenna_gain(...
                        azimuth_deg, i);
                end
                a = a .* antenna_gains;
            end

            % Normalize
            a = a / norm(a);
        end

        function R = compute_covariance_matrix(obj, signals, ...
                use_forward_backward)
            % Compute spatial covariance matrix
            %
            % Inputs:
            %   signals: (M×N) matrix, M=6 antennas, N=snapshots
            %   use_forward_backward: Use FB averaging (default true)
            %
            % Returns:
            %   R: (M×M) covariance matrix

            if nargin < 3
                use_forward_backward = true;
            end

            [M, N] = size(signals);

            % Forward covariance
            R_forward = (1/N) * (signals * signals');

            if use_forward_backward
                % Exchange matrix
                J = eye(M);
                J = J(end:-1:1, :);

                % Backward covariance
                R_backward = J * conj(R_forward) * J;

                % Average
                R = 0.5 * (R_forward + R_backward);
            else
                R = R_forward;
            end
        end

        function [eigenvalues, eigenvectors] = eigendecomposition(~, R)
            % Eigendecomposition of covariance matrix
            %
            % Inputs:
            %   R: (M×M) covariance matrix
            %
            % Returns:
            %   eigenvalues: (M×1) eigenvalues (descending order)
            %   eigenvectors: (M×M) eigenvectors

            % Use eig for Hermitian matrix
            [V, D] = eig(R);

            % Sort by magnitude (descending)
            [eigenvalues, idx] = sort(diag(D), 'descend');
            eigenvectors = V(:, idx);
        end

        function num_src = estimate_num_sources(obj, eigenvalues, method)
            % Estimate number of sources using MDL or AIC
            %
            % Inputs:
            %   eigenvalues: Eigenvalue array
            %   method: 'MDL' (default) or 'AIC'
            %
            % Returns:
            %   num_src: Estimated number of sources

            if nargin < 3
                method = 'MDL';
            end

            M = length(eigenvalues);
            N = 1000;  % Assumed snapshots (should be passed as parameter)

            if strcmpi(method, 'MDL')
                mdl = zeros(M, 1);
                for k = 0:M-1
                    noise_eigs = eigenvalues(k+1:end);
                    geo_mean = exp(mean(log(noise_eigs + 1e-12)));
                    arith_mean = mean(noise_eigs);

                    mdl(k+1) = -N * (M - k) * log(geo_mean / arith_mean) + ...
                        0.5 * k * (2*M - k) * log(N);
                end
                [~, num_src] = min(mdl);
                num_src = num_src - 1;  % Adjust for 0-indexing

            else  % AIC
                aic = zeros(M, 1);
                for k = 0:M-1
                    noise_eigs = eigenvalues(k+1:end);
                    geo_mean = exp(mean(log(noise_eigs + 1e-12)));
                    arith_mean = mean(noise_eigs);

                    aic(k+1) = -2 * N * (M - k) * log(geo_mean / arith_mean) + ...
                        2 * k * (2*M - k);
                end
                [~, num_src] = min(aic);
                num_src = num_src - 1;
            end

            fprintf('Estimated %d sources using %s\n', num_src, method);
        end

        function result = compute_music_spectrum(obj, R, ...
                azimuth_range, resolution, auto_estimate_sources)
            % Compute MUSIC angular spectrum
            %
            % Inputs:
            %   R: Covariance matrix
            %   azimuth_range: [min, max] in degrees (default [-180, 180])
            %   resolution: Angular resolution in degrees (default 1.0)
            %   auto_estimate_sources: Auto-estimate num sources (default false)
            %
            % Returns:
            %   result: struct with fields:
            %     - azimuth: Azimuth grid
            %     - spectrum: MUSIC spectrum (dB)
            %     - eigenvalues: Eigenvalues
            %     - num_sources: Number of sources used

            if nargin < 3
                azimuth_range = [-180, 180];
            end
            if nargin < 4
                resolution = 1.0;
            end
            if nargin < 5
                auto_estimate_sources = false;
            end

            % Eigendecomposition
            [eigenvalues, eigenvectors] = obj.eigendecomposition(R);

            % Auto-estimate number of sources
            num_sources = obj.num_sources;
            if auto_estimate_sources
                num_sources = obj.estimate_num_sources(eigenvalues);
            end

            % Noise subspace
            noise_subspace = eigenvectors(:, num_sources+1:end);

            % Create angular grid
            azimuth_grid = azimuth_range(1):resolution:azimuth_range(2);

            % Compute MUSIC spectrum (vectorized)
            spectrum = zeros(length(azimuth_grid), 1);

            for i = 1:length(azimuth_grid)
                az = azimuth_grid(i);

                % Steering vector
                a = obj.steering_vector(az, 0);

                % MUSIC pseudo-spectrum: P(θ) = 1 / ||En^H * a||^2
                projection = noise_subspace' * a;
                denominator = abs(projection' * projection);

                spectrum(i) = 1.0 / (denominator + 1e-12);
            end

            % Convert to dB (normalized)
            spectrum_db = 10 * log10(spectrum / max(spectrum) + 1e-12);

            % Return results
            result = struct();
            result.azimuth = azimuth_grid;
            result.spectrum = spectrum_db;
            result.eigenvalues = eigenvalues;
            result.num_sources = num_sources;
        end

        function detected_sources = find_peaks(obj, spectrum_result, ...
                num_peaks, min_peak_height, min_peak_separation)
            % Find peaks in MUSIC spectrum
            %
            % Inputs:
            %   spectrum_result: Output from compute_music_spectrum()
            %   num_peaks: Number of peaks (default: num_sources)
            %   min_peak_height: Min height in dB (default: -10)
            %   min_peak_separation: Min separation in degrees (default: 10)
            %
            % Returns:
            %   detected_sources: Array of structs with fields:
            %     - azimuth, elevation, magnitude_db, confidence, type

            if nargin < 3 || isempty(num_peaks)
                num_peaks = spectrum_result.num_sources;
            end
            if nargin < 4
                min_peak_height = -30.0;  % Further lowered threshold for weak peaks
            end
            if nargin < 5
                min_peak_separation = 30.0;  % Increased to prevent ghost peaks (180° vs 189°)
            end

            spectrum = spectrum_result.spectrum;
            azimuth = spectrum_result.azimuth;

            % Find local maxima
            resolution = mean(diff(azimuth));
            window_size = max(3, round(min_peak_separation / resolution));

            local_max = movmax(spectrum, window_size, 'Endpoints', 'fill');

            % Peaks are local maxima above threshold
            is_peak = (spectrum == local_max) & (spectrum > min_peak_height);

            peak_indices = find(is_peak);

            if isempty(peak_indices)
                warning('No peaks found in MUSIC spectrum');
                detected_sources = struct([]);
                return;
            end

            % Sort by magnitude
            peak_values = spectrum(peak_indices);
            [~, sorted_idx] = sort(peak_values, 'descend');

            % Extract top peaks
            detected_sources = struct([]);
            for i = 1:min(num_peaks, length(sorted_idx))
                idx = sorted_idx(i);
                peak_idx = peak_indices(idx);

                % Parabolic interpolation for sub-degree accuracy
                if peak_idx > 1 && peak_idx < length(spectrum)
                    y1 = spectrum(peak_idx - 1);
                    y2 = spectrum(peak_idx);
                    y3 = spectrum(peak_idx + 1);
                    delta = 0.5 * (y1 - y3) / (y1 - 2*y2 + y3 + 1e-12);
                    az_refined = azimuth(peak_idx) + delta * resolution;
                else
                    az_refined = azimuth(peak_idx);
                end

                source.azimuth = az_refined;
                source.elevation = 0.0;
                source.magnitude_db = spectrum(peak_idx);
                source.confidence = (spectrum(peak_idx) - min_peak_height) / ...
                    (0 - min_peak_height);
                source.type = 'MUSIC_peak';

                detected_sources = [detected_sources; source];

                fprintf('Detected source #%d: Az=%.2f°, Mag=%.1f dB, Conf=%.2f\n', ...
                    i, source.azimuth, source.magnitude_db, source.confidence);
            end
        end

        function [detected_sources, spectrum_result] = estimate_aoa(obj, ...
                signals, azimuth_range, resolution, use_forward_backward, ...
                auto_estimate_sources)
            % Complete MUSIC AOA estimation pipeline
            %
            % Inputs:
            %   signals: (6×N) signal data array
            %   azimuth_range: [min, max] (default [-180, 180])
            %   resolution: Angular resolution (default 0.5)
            %   use_forward_backward: Use FB averaging (default true)
            %   auto_estimate_sources: Auto-estimate sources (default false)
            %
            % Returns:
            %   detected_sources: Array of detected source structs
            %   spectrum_result: MUSIC spectrum result (for visualization)

            if nargin < 3
                azimuth_range = [-180, 180];
            end
            if nargin < 4
                resolution = 0.5;
            end
            if nargin < 5
                use_forward_backward = true;
            end
            if nargin < 6
                auto_estimate_sources = false;
            end

            % Validate input
            if size(signals, 1) ~= obj.num_elements
                error('Expected %d antennas, got %d', ...
                    obj.num_elements, size(signals, 1));
            end

            fprintf('Starting AOA estimation: %d snapshots, resolution=%.1f°\n', ...
                size(signals, 2), resolution);

            % Compute covariance matrix
            R = obj.compute_covariance_matrix(signals, use_forward_backward);

            % Compute MUSIC spectrum
            spectrum_result = obj.compute_music_spectrum(R, ...
                azimuth_range, resolution, auto_estimate_sources);

            % Find peaks
            detected_sources = obj.find_peaks(spectrum_result);

            fprintf('AOA estimation complete: %d sources detected\n', ...
                length(detected_sources));
        end
    end
end
