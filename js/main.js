import * as param from "./const.js"
import paramIter from "./paramIter.js";
import { vk, rk } from "./express.js";

// let i = 1;
const integrXLSX = (thet, h, t1, t2, P0, step = 1) => {
    const body = document.querySelector('body')
    let dt = param.step;
    let t = param.t0;
    let el = [param.m0, param.vx0, param.vy0, param.x0, param.y0, param.v0, param.h0, param.r0, param.thet * 180 / math.pi, param.THETc * 180 / math.pi, param.THET * 180 / math.pi, param.alfa * 180 / math.pi, param.fi * 180 / math.pi];

    let t_prev = t;
    let el_prev = el.slice();
    let boolk = false;
    body.insertAdjacentHTML('beforeend', `<div>t1 = ${t1}</div>`);
    body.insertAdjacentHTML('beforeend', `<div>t2 = ${t2}</div>`);
    body.insertAdjacentHTML('beforeend', `<div>thet_torch = ${thet[0]}</div>`);
    body.insertAdjacentHTML('beforeend', `<div>thet2 = ${thet[1]}</div>`);
    body.insertAdjacentHTML('beforeend', `<div> t, [-]   |m, [м]     |vx, [м/c]     |vy [м/c]     |x, [м]     |y, [м]     |v, [м/с]     |h, [м]    |r, [м]     |ϑ, [град]     |θс, [град]     |θ, [град]     |α, [град]     |φ, [град]     </div>`);
    body.insertAdjacentHTML('beforeend', `<div>${t.toFixed(5)}|${el[0].toFixed(5)}|${el[1].toFixed(5)}|${el[2].toFixed(5)}|${el[3].toFixed(5)}|${el[4].toFixed(5)}|${el[5].toFixed(5)}|${el[6].toFixed(5)}|${el[7].toFixed(5)}|${el[8].toFixed(5)}|${el[9].toFixed(5)}|${el[10].toFixed(5)}|${el[11].toFixed(5)}|${el[12].toFixed(5)}</div>`);
    while (math.abs(el[5] - vk(h)) > param.eps_v) {
        if (t + dt > t1 && t < t1) {
            dt = t1 - t;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0, boolk, true);
            // body.insertAdjacentHTML('beforeend', `<div>${t.toFixed(5)}|${el[0].toFixed(5)}|${el[1].toFixed(5)}|${el[2].toFixed(5)}|${el[3].toFixed(5)}|${el[4].toFixed(5)}|${el[5].toFixed(5)}|${el[6].toFixed(5)}|${el[7].toFixed(5)}|${el[8].toFixed(5)}|${el[9].toFixed(5)}|${el[10].toFixed(5)}|${el[11].toFixed(5)}|${el[12].toFixed(5)}</div>`);
            t += dt;
            boolk = true;
            dt = param.step - dt;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0, boolk, true);
            body.insertAdjacentHTML('beforeend', `<div>${t.toFixed(5)}|${el[0].toFixed(5)}|${el[1].toFixed(5)}|${el[2].toFixed(5)}|${el[3].toFixed(5)}|${el[4].toFixed(5)}|${el[5].toFixed(5)}|${el[6].toFixed(5)}|${el[7].toFixed(5)}|${el[8].toFixed(5)}|${el[9].toFixed(5)}|${el[10].toFixed(5)}|${el[11].toFixed(5)}|${el[12].toFixed(5)}</div>`);
            t += dt;
            dt = param.step;
        }
        if (t + dt > t2 && t < t2) {
            dt = t2 - t;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0, boolk, true);
            // body.insertAdjacentHTML('beforeend', `<div>${t.toFixed(5)}|${el[0].toFixed(5)}|${el[1].toFixed(5)}|${el[2].toFixed(5)}|${el[3].toFixed(5)}|${el[4].toFixed(5)}|${el[5].toFixed(5)}|${el[6].toFixed(5)}|${el[7].toFixed(5)}|${el[8].toFixed(5)}|${el[9].toFixed(5)}|${el[10].toFixed(5)}|${el[11].toFixed(5)}|${el[12].toFixed(5)}</div>`);
            t += dt;
            boolk = false;
            dt = param.step - dt;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0, boolk, true);
            body.insertAdjacentHTML('beforeend', `<div>${t.toFixed(5)}|${el[0].toFixed(5)}|${el[1].toFixed(5)}|${el[2].toFixed(5)}|${el[3].toFixed(5)}|${el[4].toFixed(5)}|${el[5].toFixed(5)}|${el[6].toFixed(5)}|${el[7].toFixed(4)}|${el[8].toFixed(5)}|${el[9].toFixed(5)}|${el[10].toFixed(5)}|${el[11].toFixed(5)}|${el[12].toFixed(5)}</div>`);
            t += dt;
            dt = param.step
        }
        el_prev = el.slice();
        t_prev = t;
        paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0, boolk, true);
        if ((el[5] < vk(h)) && (+t.toFixed(1) % step === 0 && +t.toFixed(1) !== 0)) {
            body.insertAdjacentHTML('beforeend', `<div>${t.toFixed(5)}|${el[0].toFixed(5)}|${el[1].toFixed(5)}|${el[2].toFixed(5)}|${el[3].toFixed(5)}|${el[4].toFixed(5)}|${el[5].toFixed(5)}|${el[6].toFixed(5)}|${el[7].toFixed(5)}|${el[8].toFixed(5)}|${el[9].toFixed(5)}|${el[10].toFixed(5)}|${el[11].toFixed(5)}|${el[12].toFixed(5)}</div>`);
        }
        if (math.abs(el[5] - vk(h)) < param.eps_v) {
            body.insertAdjacentHTML('beforeend', `<div>${t.toFixed(5)}|${el[0].toFixed(5)}|${el[1].toFixed(5)}|${el[2].toFixed(5)}|${el[3].toFixed(5)}|${el[4].toFixed(5)}|${el[5].toFixed(5)}|${el[6].toFixed(5)}|${el[7].toFixed(5)}|${el[8].toFixed(5)}|${el[9].toFixed(5)}|${el[10].toFixed(5)}|${el[11].toFixed(5)}|${el[12].toFixed(5)}</div>`);
        }
        t += dt;


        if (el[5] > vk(h)) {
            el = el_prev.slice();
            t = t_prev;
            el_prev = el.slice();
            t_prev = t;
            dt = dt / 10;
        }
    }
    return el
}
const integr = (thet, h, t1, t2, P0) => {
    let dt = param.step;
    let t = param.t0;
    let el = [param.m0, param.vx0, param.vy0, param.x0, param.y0, param.v0, param.h0, param.r0, param.thet, param.THETc, param.THET, param.alfa, param.fi];

    let t_prev = t;
    let el_prev = el.slice();
    let boolk = false;
    while (math.abs(el[5] - vk(h)) > param.eps_v) {
        if (t + dt > t1 && t < t1) {
            dt = t1 - t;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0, boolk);
            t += dt;
            boolk = true;
            dt = param.step - dt;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0, boolk);
            t += dt;
            dt = param.step;
        }
        if (t + dt > t2 && t < t2) {
            dt = t2 - t;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0, boolk);
            t += dt;
            boolk = false;
            dt = param.step - dt;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0, boolk);
            t += dt;
            dt = param.step
        }
        el_prev = el.slice();
        t_prev = t;
        paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0, boolk);
        t += dt;

        if (el[5] > vk(h)) {
            el = el_prev.slice();
            t = t_prev;
            el_prev = el.slice();
            t_prev = t;
            dt = dt / 10;
        }
    }
    return el
}
const optimus = (thet_arr, t, P0, h) => {
    let thet = thet_arr.slice();
    const thet_step = [math.pow(10, -5), math.pow(10, -8)];
    let el1, dF;
    function thet_1_2_torch(thet) {
        let el, thetk, delta_next_r_1, delta_next_r_2, delta_next_thet_1, delta_next_thet_2;
        let delta_prev_r_1, delta_prev_r_2, delta_prev_thet_1, delta_prev_thet_2;
        thetk = thet.slice();
        thetk[0] += thet_step[0];
        el = integr(thetk, h, t[0], t[1], P0);
        delta_next_r_1 = el[7] - rk(h);
        delta_next_thet_1 = el[10];

        thetk[0] -= 2 * thet_step[0];
        el = integr(thetk, h, t[0], t[1], P0);
        delta_prev_r_1 = el[7] - rk(h);
        delta_prev_thet_1 = el[10];

        thetk = thet.slice();
        thetk[1] += thet_step[1];
        el = integr(thetk, h, t[0], t[1], P0);
        delta_next_r_2 = el[7] - rk(h);
        delta_next_thet_2 = el[10];

        thetk[1] -= 2 * thet_step[1];
        el = integr(thetk, h, t[0], t[1], P0);
        delta_prev_r_2 = el[7] - rk(h);
        delta_prev_thet_2 = el[10];

        const dr_1 = (delta_next_r_1 - delta_prev_r_1) / (2 * thet_step[0]);
        const dthet_1 = (delta_next_thet_1 - delta_prev_thet_1) / (2 * thet_step[0]);
        const dr_2 = (delta_next_r_2 - delta_prev_r_2) / (2 * thet_step[1]);
        const dthet_2 = (delta_next_thet_2 - delta_prev_thet_2) / (2 * thet_step[1]);
        return [
            [dr_1, dr_2],
            [dthet_1, dthet_2]
        ]
    }
    do {
        el1 = integr(thet, h, t[0], t[1], P0);
        dF = [-el1[7] + rk(h), -el1[10]];
        if (math.sqrt(math.pow((dF[0]) / param.eps_r, 2) + math.pow((dF[1]) / param.eps_thet, 2)) < 1) {
            break
        }
        const J = math.matrix(thet_1_2_torch(thet));
        const J_inv = math.inv(J);
        const U = math.multiply(J_inv, dF);
        thet = math.add(thet, U)._data;
    } while (true)
    return [
        thet,
        el1[0]
    ]
}

function levenberg(P, h, t) {
    function dm_t1(t0_arr) {
        let delta_next, delta_prev, tk_arr;
        tk_arr = t0_arr.slice();
        tk_arr[0] += dt_step;
        delta_next = optimus(thet_id, tk_arr, P, h)[1];
        tk_arr[0] -= 2 * dt_step;
        delta_prev = optimus(thet_id, tk_arr, P, h)[1];
        return (delta_next - delta_prev) / (2 * dt_step)
    }

    function dm_t2(t0_arr) {
        let delta_next, delta_prev, tk_arr;
        tk_arr = t0_arr.slice();
        tk_arr[1] += dt_step;
        delta_next = optimus(thet_id, tk_arr, P, h)[1];
        tk_arr[1] -= 2 * dt_step;
        delta_prev = optimus(thet_id, tk_arr, P, h)[1];
        return (delta_next - delta_prev) / (2 * dt_step)
    }

    function ddm_t1(t0_arr) {
        let z1, z2, f_res, dm1_1, dm1_2, dm2_1, dm2_2, tk_arr;
        tk_arr = t0_arr.slice();
        dm1_2 = optimus(thet_id, tk_arr, P, h)[1];
        dm2_1 = dm1_2;
        tk_arr[0] += 2 * dt_step_double;
        dm2_2 = optimus(thet_id, tk_arr, P, h)[1];
        tk_arr[0] -= 4 * dt_step_double;
        dm1_1 = optimus(thet_id, tk_arr, P, h)[1];
        z1 = (dm1_1 - dm1_2) / (2 * dt_step_double);
        z2 = (dm2_1 - dm2_2) / (2 * dt_step_double);
        f_res = (z1 - z2) / (2 * dt_step_double);
        return f_res
    }

    function ddm_t2(t0_arr) {
        let z1, z2, f_res, dm1_1, dm1_2, dm2_1, dm2_2, tk_arr;
        tk_arr = t0_arr.slice();
        dm1_2 = optimus(thet_id, tk_arr, P, h)[1];
        dm2_1 = dm1_2;
        tk_arr[1] += 2 * dt_step_double;
        dm2_2 = optimus(thet_id, tk_arr, P, h)[1];
        tk_arr[1] -= 4 * dt_step_double;
        dm1_1 = optimus(thet_id, tk_arr, P, h)[1];
        z1 = (dm1_1 - dm1_2) / (2 * dt_step_double);
        z2 = (dm2_1 - dm2_2) / (2 * dt_step_double);
        f_res = (z1 - z2) / (2 * dt_step_double);
        return f_res
    }
    function ddm_t1_t2(t0_arr) {
        let z1, z2, f_res, dm1_1, dm1_2, dm2_1, dm2_2, tk_arr;
        tk_arr = t0_arr.slice();
        tk_arr[0] += dt_step;
        tk_arr[1] += dt_step_double;
        dm1_1 = optimus(thet_id, tk_arr, P, h)[1];
        tk_arr[0] -= 2 * dt_step;
        dm2_2 = optimus(thet_id, tk_arr, P, h)[1];
        tk_arr[0] += 2 * dt_step;
        tk_arr[1] -= 2 * dt_step_double;
        dm2_1 = optimus(thet_id, tk_arr, P, h)[1];
        tk_arr[0] -= 2 * dt_step;
        dm1_2 = optimus(thet_id, tk_arr, P, h)[1];
        z1 = (dm1_1 - dm2_2) / (2 * dt_step);
        z2 = (dm2_1 - dm1_2) / (2 * dt_step);
        f_res = (z1 - z2) / (2 * dt_step_double);
        return f_res
    }

    let t_id = t.slice();
    let dt_step = math.pow(10, -2);
    let dt_step_double = 2 * math.pow(10, -2);

    const I = math.matrix(
        [
            [1, 0],
            [0, 1]
        ]
    )
    let C1 = math.pow(2, -1);
    let C2 = math.pow(2, 1);
    let alfa = math.pow(10, -4);
    let alfa_iter = alfa;
    let t_iter = t_id.slice();
    let t_prev = t_iter.slice();

    let gradient, gradient_module
    let m_iter, m_iter_prev
    let thet_id = optimus([param.thet_torch, param.thet_2], t_iter, P, h)[0];
    do {
        gradient = math.matrix([dm_t1(t_iter), dm_t2(t_iter)]);
        gradient_module = math.sqrt(math.pow(gradient._data[0], 2) + math.pow(gradient._data[1], 2));
        m_iter = optimus(thet_id, t_iter, P, h)[1];
        console.log("Время: ", t_iter);
        console.log("Градиент: ", gradient_module);
        console.log("Mass", m_iter);
        console.log("alfa", alfa_iter);
        if (gradient_module < param.eps_extr) {
            break
        }
        // Гессиан в цикле пока модуль больше предыдушего значения модуля
        // gradient_prev = math.clone(gradient);
        m_iter_prev = m_iter;
        t_prev = t_iter.slice();
        let gessian = math.matrix(
            [
                [ddm_t1(t_iter), ddm_t1_t2(t_iter)],
                [ddm_t1_t2(t_iter), ddm_t2(t_iter)]
            ]
        )
        do {
            let first = math.add(gessian, math.multiply(alfa_iter, I)); // (Hi + alfai * I)
            let second = math.inv(first); // // (Hi + alfai * I)^-1
            let third = math.multiply(second, gradient); // (Hi + alfai * I)^-1 x gradient
            // third = math.multiply(third, -1);
            t_iter = math.add(t_iter, third)._data;
            // gradient = math.matrix([dm_t1(t_iter), dm_t2(t_iter)]);
            // gradient_module = math.sqrt(math.pow(gradient._data[0], 2) + math.pow(gradient._data[1], 2));
            m_iter = optimus(thet_id, t_iter, param.P2, param.h_isl_2_2)[1];

            if (m_iter > m_iter_prev) {
                alfa_iter *= C1;
                break;
            }
            alfa_iter *= C2;
            t_iter = t_prev.slice();
        } while (true);
    } while (true);
    return t_iter
}
function init() {
    // let t_id = [360, 700]
    // levenberg(param.P1, param.h_isl_1_2, t_id)

    // МОЕ
    let thet_id = [param.thet_torch, param.thet_2]
    let stepXLSX = 50;
    let t1 = [481.7201771642, 749.0503961927]
    // let thet1 = optimus(thet_id, t1, param.P1, param.h_isl_1_1)[0]
    // let traekt1 = integrXLSX(thet1, param.h_isl_1_1, t1[0], t1[1], param.P1)

    let t2 = [524.4050296731, 978.4411379483]
    // let thet2 = optimus(thet_id, t2, param.P1, param.h_isl_1_2)[0]
    // let traekt2 = integrXLSX(thet2, param.h_isl_1_2, t2[0], t2[1], param.P1)

    let t3 = [396.6877195853, 792.2535392132]
    // let thet3 = optimus(thet_id, t3, param.P2, param.h_isl_2_1)[0]
    // let traekt3 = integrXLSX(thet3, param.h_isl_2_1, t3[0], t3[1], param.P2)

    let t4 = [460.7020301335, 1426.1057798270]
    // let thet4 = optimus(thet_id, t4, param.P2, param.h_isl_2_2)[0]
    // let traekt4 = integrXLSX(thet4, param.h_isl_2_2, t4[0], t4[1], param.P2)
    ///////////////
    function t1k(t1, h, P, sigma, chad) {
        const body = document.querySelector('body')
        let tk = t1[0] - sigma;
        let dt1 = (t1[1] - t1[0]) / chad;
        let mk = 0;
        body.insertAdjacentHTML('beforeend', `<div>t1, [c] | mk, [кг]</div>`);
        while (tk <= t1[1]) {
            mk = optimus(thet_id, [tk, t1[1]], P, h)[1]
            body.insertAdjacentHTML('beforeend', `<div>${tk} | ${mk}</div>`);
            tk += dt1
        }
    }
    function t2k(t1, h, P, chad) {
        const body = document.querySelector('body')
        let tk = 1500;
        let dt1 = (tk - t1[0]) / chad;
        let mk = 0;
        body.insertAdjacentHTML('beforeend', `<div>t1, [c] | mk, [кг]</div>`);
        do {
            mk = optimus(thet_id, [t1[0], tk], P, h)[1]
            console.log(mk);
            body.insertAdjacentHTML('beforeend', `<div>${tk} | ${mk}</div>`);
            tk -= dt1
        } while (tk >= t1[0])
    }
    // t1k(t1, param.h_isl_1_1, param.P1,-50, 15)
    // t1k(t2, param.h_isl_1_2, param.P1,-50, 15)
    // t1k(t3, param.h_isl_2_1, param.P2, -50, 15)
    // t1k(t4, param.h_isl_2_2, param.P2, -100, 15)

    // let m1 = optimus(thet_id, t3, param.P2, param.h_isl_1_2)[1] - param.m_const
    // console.log(m1);
    function orbita_param(Vx0, Vy0, x0, y0) {
        // debugger
        let mu = 4.903 * 1e12;
        let z0 = 0;
        let Vz0 = 0;
        let r0 = math.sqrt(math.pow(x0, 2) + math.pow(y0 + param.rM, 2) + math.pow(z0, 2));
        let V0 = math.sqrt(math.pow(Vx0, 2) + math.pow(Vy0, 2) + math.pow(Vz0, 2));
        let k0 = r0 * math.pow(V0, 2) / mu;
        let tetta0 = math.asin((x0 * Vx0 + (y0 + param.rM) * Vy0 + z0 * Vz0) / (r0 * V0)); // угол наклонения вектора скорости
        let sintetta0 = (x0 * Vx0 + (y0 + param.rM) * Vy0 + z0 * Vz0) / (r0 * V0);
        let costetta0 = (math.sqrt(math.pow(r0, 2) * math.pow(V0, 2) - math.pow(x0 * Vx0 + (y0 + param.rM) * Vy0 + z0 * Vz0, 2))) / (r0 * V0);
        if (tetta0 < -math.pi / 2 && tetta0 > math.pi / 2) {
            console.log("error");
        }
        let a = 0;   // большая полуось
        let e = 0;   // эксцентриситет
        let omg = 0; // долгота восходящего узла
        let tau = 0; // время прохождения через перигей
        let w = 0;   // аргумент перицентра
        let TA = 0;  // истинная аномалия

        let u0 = 0;
        a = r0 / (2 - k0);
        e = math.sqrt(math.pow(1 - k0, 2) + k0 * (2 - k0) * math.pow(sintetta0, 2));
        
        let C1 = y0 * Vz0 - z0 * Vy0;
        let C2 = z0 * Vx0 - x0 * Vz0;
        let C3 = x0 * Vy0 - y0 * Vx0;
        let C = math.sqrt(math.pow(C1, 2) + math.pow(C2, 2) + math.pow(C3, 2));

        let i = math.acos(C3 / C);

        if (i < 0 && i > math.pi) {
            console.log("error");
        }
        if (math.abs(i - 0) < 1e-6) {
            omg = 0;
            w = 0;
            if (x0 < 0 || y0 < 0 || z0 < 0) {
                u0 = math.acos((x0 * math.cos(omg) + y0 * math.sin(omg)) / r0);
            }
            else {
                u0 = math.acos((x0 * cos(omg + math.pi) + y0 * math.sin(omg + math.pi)) / r0);
            }
            TA = u0;
        }
        else {
            if (x0 < 0 || y0 < 0 || z0 < 0 || math.abs(x0 - 0) < 1e-6 || math.abs(y0 - 0) < 1e-6 || math.abs(z0 - 0) < 1e-6) {
                omg = math.atan(-C1 / C2);
            }
            else {
                omg = math.atan(-C1 / C2) + math.pi;
            }
            if (e < 1e-6) {
                if (x0 < 0 || y0 < 0 || z0 < 0 || math.abs(x0 - 0) < 1e-6 || math.abs(y0 - 0) < 1e-6 || math.abs(z0 - 0) < 1e-6) {
                    u0 = math.acos((x0 * math.cos(omg) + y0 * math.sin(omg)) / r0);
                }
                else {
                    u0 = math.acos((x0 * math.cos(omg + math.pi) + y0 * math.sin(omg + math.pi)) / r0);
                }
                TA = u0;
                w = u0 - TA;
            }
            else {
                let sinTA = (k0 * sintetta0 * costetta0) / e;
                let cosTA = ((k0 * math.pow(costetta0, 2)) - 1) / e;
                let k8 = math.pow(sinTA, 2) + math.pow(cosTA, 2);
                TA = math.asin((k0 * sintetta0 * costetta0) / e);
                let TA1 = math.acos((k0 * math.pow(costetta0, 2) - 1) / e);
                let u0 = 0;
                if (x0 < 0 || y0 < 0 || z0 < 0 || math.abs(x0 - 0) < 1e-6 || math.abs(y0 - 0) < 1e-6 || math.abs(z0 - 0) < 1e-6) {
                    let op = ((x0 * math.cos(omg) + y0 * math.sin(omg)) / r0);
                    if (math.abs(op - 1) < 1e-6) {
                        op = 1;
                        u0 = math.acos(op);
                    }
                    else {
                        u0 = math.acos((x0 * math.cos(omg) + y0 * math.sin(omg)) / r0);
                    }
                }
                else {
                    u0 = math.acos((x0 * math.cos(omg + math.pi) + y0 * math.sin(omg + math.pi)) / r0);
                }
                w = u0 - TA;
            }
        }
        i = i * 180 / math.pi;
        if (omg * 180 / math.pi < 0) {
            omg = 360 + omg * 180 / math.pi;
        }
        else {
            omg *= 180 / math.pi;
        }
        w *= 180 / math.pi;
        if (TA * 180 / math.pi < 0) {
            TA = 360 + TA * 180 / math.pi;
        }
        else {
            TA *= 180 / math.pi;
        }
        let ra = a * (1 + e);
        let rp = a * (1 - e);
        const body = document.querySelector('body')
        body.insertAdjacentHTML('beforeend', `<div>${rp}</div>`)
        body.insertAdjacentHTML('beforeend', `<div>${ra}</div>`)
        body.insertAdjacentHTML('beforeend', `<div>${a}</div>`)
        body.insertAdjacentHTML('beforeend', `<div>${e}</div>`)
        return [
            rp,
            ra,
            a,
            e
            // omg,
            // tau,
            // w,
            // TA,
            // i,
        ];
    }
    // const orbita = orbita_param(1166.54262,	109.72362	,189161.2462,	93920.3863)
    // const orbita = orbita_param(1333.28224,	69.66524,	236460.9296,	106633.0777)
    // const orbita = orbita_param(1134.46076,	239.01589,	150197.7036,	101373.1737)
    const orbita = orbita_param(1441.7898,	169.354	,224120.8207	,127385.608)
    console.log(orbita);
    // ДАНИК
    // let P1 = 8.6 * math.pow(10,3);
    // let h1 = 161 * math.pow(10,3);
    // let h2 = 205 * math.pow(10,3);
    // let t1 = [400.1221, 700.49494]
    // let thet1 = optimus(thet_id, t1, P1, h1)[0];
    // console.log(thet1);
    // let traekt1 = integrXLSX(thet1, h1, t1[0], t1[1], P1);


}
init()
